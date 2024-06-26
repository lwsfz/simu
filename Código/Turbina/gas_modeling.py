"""
Se emplea la librería PyroMat para modelar los productos como mezcla de gases ideales no perfectos y determinar
sus propiedades a partir de las que se proporcionen.
"""
import math

import pyromat as pm
from math import fabs, sqrt

pm.config['unit_pressure'] = 'Pa'


class mixpm:
    """ Esta clase define un objeto que modela los productos de una combustión completa. Las propiedades termodinámicas
     de los componentes de la mezcla se van a calcular mediante el empleo de la librería PyroMat"""

    def __init__(self, mode: str):
        """ :param mode: Cadena de caracteres que indica el modo de funcionamiento según si se emplea, o no, modelo de
                         gas ideal en todos los componentes de los productos. Debe ser "ig" o "mp". """

        self.xmi = []  # xmi son las fracciones másicas de los componentes
        self.xni = []  # xni son las fracciones molares de los componentes
        self.mode = mode
        self._cpmix = 0.0  # Calor específico a presión constante de la mezcla
        self._cvmix = 0.0  # Calor específico a volumen constante de la mezcla
        self._gammamix = 0.0  # Coeficiente adiabático de la mezcla
        self.__Mwtotal = 0.0  # Molar weight de la mezcla
        self.mw_elem = (44.010, 18.0153, 28.0134, 31.999)
        self.sumprops_Tpd = 0.0  # Almacenaje de la parte o la totalidad de la propiedad que se calcula
        self.x = 0.0  # Átomos de carbono en cada átomo de hidrocarburo en los reactivos
        self.y = 0.0  # Átomos de hidrógeno en cada átomo de hidrocarburo en los reactivos
        self.N = 0.0  # Exceso de aire que se considera en el ajuste estequiométrico de la combustión completa
        self.mu_elem = []
        self.mu_mix = 0.0
        self.rel_error = 1E-12  # Máximo error relativo permitido.

        ''' En las líneas que siguen se almacenan en una lista atributo de la clase 'mixpm' identificadores de 
            los componentes como objetos de la librería PyroMat según el modelo especificado. '''

        modes = ("ig", "mp", )
        componentes = ("CO2", "H2O", "N2", "O2", )
        self.components = []
        if mode in modes:
            for comp in componentes:
                if comp == "H20":
                    self.components.append(pm.get(f'{mode}.{comp}'))
                else:
                    self.components.append(pm.get(f'ig.{comp}'))
            self.components = tuple(self.components)

    def setmix(self, x: float, y: float, N: float):
        """ Se determinan las propiedades del objeto que contienen las fracciones molares (xni) y másicas (xmi) de los
        productos de una combustión completa del hidrocarburo indicado.
                :param x: Átomos de carbono en cada átomo de hidrocarburo en los reactivos.
                :param y: Átomos de hidrógeno en cada átomo de hidrocarburo en los reactivos.
                :param N: Exceso de aire que se considera en el ajuste estequiométrico.
                        :return: No se devuelve nada. """
        self.x, self.y, self.N = x, y, N
        # den es la masa total por kmol de combustible en combustión completa
        den = (x * self.components[0].mw()) + (y * self.components[1].mw() / 2) + \
              (self.components[2].mw() * N * (x + y / 4) * 79.0 / 21.0) + (
                      self.components[3].mw() * (N - 1) * (x + y / 4))
        # En las líneas que siguen el factor 'f' en 'f/den' es la masa del componente por kmol keroseno
        self.xmi.append(self.components[0].mw() * x / den)
        self.xmi.append(self.components[1].mw() * y / (2 * den))
        self.xmi.append(self.components[2].mw() * N * (x + y / 4) * (79.0 / 21.0) / den)
        self.xmi.append(self.components[3].mw() * (N - 1) * (x + y / 4) / den)
        self.xmi = tuple(self.xmi)
        tmwden = 0
        for i, xm in enumerate(self.xmi):
            tmwden += xm / self.components[i].mw()
        self.__Mwtotal = 1 / tmwden
        for i, xm in enumerate(self.xmi):
            self.xni.append(xm * self.__Mwtotal / self.components[i].mw())
        return

    def get_coeffs(self, T: float, p= float):
        """ Se calculan los calores específicos y el coeficiente adiabático correspondientes con el estado de la mezcla
            que se determine.
                :param T: Temperatura de los productos (K).
                :param p: Presión total de los productos (Pa).
                        :return: En el mismo orden: Cp, Cv y gamma. Las unidades de los calores
                                específicos son kJ/(kg*K). """

        self._cpmix = 0.0
        self._cvmix = 0.0

        if self.mode == "mp":
            for i, xm in enumerate(self.xmi):  # i = 0, 1, 2, 3
                if i == 1:  # cp y cv se reciben como kJ/(kg*K)
                    self._cpmix = self._cpmix + xm * self.components[i].cp(T=T, p=p * xm)
                    self._cvmix = self._cvmix + xm * self.components[i].cv(T=T, p=p * xm)
                else:
                    self._cpmix = self._cpmix + xm * self.components[i].cp(T)
                    self._cvmix = self._cvmix + xm * self.components[i].cv(T)
        else:
            for i, xm in enumerate(self.xmi):
                self._cpmix = self._cpmix + xm * self.components[i].cp(T)
                self._cvmix = self._cvmix + xm * self.components[i].cv(T)

        self._gammamix = self._cpmix / self._cvmix
        return [float(self._cpmix), float(self._cvmix), float(self._gammamix)]

    def get_a(self, T: float, p= 101_325, extra=False):

        C_p, C_v, gamma = self.get_coeffs(T, p)
        theta_T = 3056 / T
        a = sqrt((C_p-C_v)*1000*T*(1+(
                 (gamma-1)/(1+((gamma-1) * ((theta_T**2)*((math.e**theta_T)/(((math.e**theta_T)-1)**2))))))))

        if not extra:
            return a
        else:
            return a, C_p, C_v, gamma

    def get_props_by_Tpd(self, known_props: dict, req_prop: str):
        """ Con este método, habiendo determinado el estado de la mezcla empleando una o dos propiedades de entre 'T',
        'p' o 'd' (Temperatura, presión o densidad), se podrá obtener cualquiera que se desconozca junto con 'h' y 's'.
                :param known_props: Un diccionario que debe contener, como clave, el carácter que represente la/s
                                   propiedad/es conocidas y, como valor, el valor de la misma.
                :param req_prop: El caracter que identifique la propiedad que se requiere.
                        :return: El valor de la propiedad que se requiere con las unidades de Pyromat. """

        frac_prop, self.sumprops_Tpd, code_get_sumando, xa, xb, xc = 0.0, 0.0, '', [], [], []
        k_0 = v_0 = k_1 = v_1 = None
        xj = [xa, xb]
        if self.mode == "ig":
            if req_prop == 'h':
                if len(known_props) == 2:
                    if list(known_props.keys())[1] == 'T':
                        k_0, v_0 = list(known_props.keys())[1], list(known_props.values())[1]
        if v_0 is None:
            k_0, v_0 = list(known_props.keys())[0], list(known_props.values())[0]
        var_list_gr = [req_prop, k_0]
        if not ('T' in known_props and req_prop == 'h' and self.mode == "ig") or self.mode == "mp":
            k_1 = list(known_props.keys())[1]
            v_1 = list(known_props.values())[1]
            var_list_gr.append(k_1)
            xj.append(xc)

        for j in range(len(xj)):
            for i in range(len(self.xmi)):
                exp, xji = 1, 1.0  # Establecido así considerando 'T'
                if var_list_gr[j] == 'p':
                    exp, xji = 1, self.xni[i]  # Ley de Dalton
                elif var_list_gr[j] == 'd':
                    exp, xji = 1, self.xmi[i]  # Masa del componente en el mismo volumen
                elif var_list_gr[j] in ['h', 's']:
                    exp, xji = -1, self.xmi[i]  # Solo interesa j=0, por lo que interpretar únicamente con exp=1
                if j == 0:
                    if var_list_gr[j] in ['p', 'd']:
                        exp = 0  # Suma de las presiones parciales y de las densidades parciales
                    elif var_list_gr[j] in ['h', 's']:
                        exp = exp * (-1)  # De esta manera resultan kg de productos en denominador y la suma es el total
                xj[j].append(xji ** exp)
        for i in range(len(self.xmi)):  # i = 0, 1, 2, 3
            if v_1 is None:
                code_get_sumando = f'self.sumprops_Tpd = xa[i]*float(self.components[i].{req_prop}({k_0}=' \
                                   f'v_0*xb[i]));'
            else:
                code_get_sumando = f'self.sumprops_Tpd = xa[i]*float(self.components[i].{req_prop}({k_0}=' \
                                   f'v_0*xb[i], {k_1}=v_1 * xc[i]));'
            exec(code_get_sumando)
            frac_prop += self.sumprops_Tpd
            if req_prop == 'T':
                break

        return frac_prop

    def get_props_with_hs(self, known_props: dict, init_guess: dict):
        """ Con este método y mediante el método numérico régula falsi, se posibilita usar como propiedad conocida la
        entropía o la entalpía de la mezcla de gases, sin embargo, no será posible establecer la otra como output
        requerido, por no poder ser entrada del método 'get_props_with_Tpd'.
                :param known_props: Diccionario que a cada clave, caracter que identifica una propiedad que se conoce,
                                   le asigna el valor de la misma.
                :param init_guess: Diccionario que indica la clave y la estimación inicial de la propiedad a determinar.
                        :return: El valor de la propiedad calculada. """

        k_0, v_0, k_1, v_1 = list(known_props.keys())[0], list(known_props.values())[0], '', 0.0
        k_a, k_b, v_a, v_b = '', '', 0.0, 0.0
        iter_counter = 0

        if len(known_props) == 2:
            k_1 = list(known_props.keys())[1]
            v_1 = list(known_props.values())[1]

        if k_0 in ['h', 's']:
            k_a, v_a = k_0, v_0
            if len(known_props) == 2:
                k_b, v_b = k_1, v_1
        else:
            k_a, v_a, k_b, v_b = k_1, v_1, k_0, v_0

        fvg = fabs(2 * v_a * self.rel_error)
        k_g, v_g, vp_g, bolz = list(init_guess.keys())[0], list(init_guess.values())[0], 0.0, None
        v_g1, v_g2, start, f1, f2, fpu = v_g*0.999, v_g*1.001, True, 0.0, 0.0, 0.0

        # Este bloque while es para que se verifique la condición del tma. de Bolzano
        while bolz is None or bolz > 0:
            iter_counter += 1
            if iter_counter > 1000:
                raise ValueError('Se ha alcanzado el número máximo de iteraciones.')
            if 'T' in init_guess and k_a == 'h' and self.mode == "ig":
                if bolz is None:
                    f1 = self.get_props_by_Tpd({k_g: v_g1}, k_a) - v_a
                    f2 = self.get_props_by_Tpd({k_g: v_g2}, k_a) - v_a
                else:
                    fp = (f2 - f1) / (v_g2 - v_g1)
                    fpu = fp/fabs(fp)
                    if f1*fpu > 0:
                        v_g2 = v_g1
                        f2 = f1
                        v_g1 = v_g1 * 0.99
                        f1 = self.get_props_by_Tpd({k_g: v_g1}, k_a) - v_a
                    elif f2*fpu < 0:
                        v_g1 = v_g2
                        f1 = f2
                        v_g2 = v_g2 * 1.01
                        f2 = self.get_props_by_Tpd({k_g: v_g2}, k_a) - v_a
            else:
                if bolz is None:
                    f1 = self.get_props_by_Tpd({k_g: v_g1, k_b: v_b}, k_a) - v_a
                    f2 = self.get_props_by_Tpd({k_g: v_g2, k_b: v_b}, k_a) - v_a
                else:
                    fp = (f2 - f1) / (v_g2 - v_g1)
                    fpu = fp/fabs(fp)
                    if f1*fpu > 0:
                        v_g2 = v_g1
                        f2 = f1
                        v_g1 = v_g1 * 0.99
                        f1 = self.get_props_by_Tpd({k_g: v_g1, k_b: v_b}, k_a) - v_a
                    elif f2*fpu < 0:
                        v_g1 = v_g2
                        f1 = f2
                        v_g2 = v_g2 * 1.01
                        f2 = self.get_props_by_Tpd({k_g: v_g2, k_b: v_b}, k_a) - v_a
            bolz = f1 * f2

        c = float(v_g2 - (f2 * (v_g2 - v_g1) / (f2 - f1)))
        if len(known_props) == 2:
            fc = self.get_props_by_Tpd({k_g: c, k_b: v_b}, k_a) - v_a
        else:
            fc = self.get_props_by_Tpd({k_g: c}, k_a) - v_a

        iter_counter = 0

        while fvg > self.rel_error:
            iter_counter += 1
            if iter_counter > 1000:
                raise ValueError('Se ha alcanzado el número máximo de iteraciones.')
            if not start:
                if 'T' in init_guess and k_a == 'h' and self.mode == "ig":
                    c = float(v_g2 - (f2 * (v_g2 - v_g1) / (f2 - f1)))
                    fc = self.get_props_by_Tpd({k_g: c}, k_a) - v_a
                else:
                    c = float(v_g2 - (f2 * (v_g2 - v_g1) / (f2 - f1)))
                    fc = self.get_props_by_Tpd({k_g: c, k_b: v_b}, k_a) - v_a
            else:
                start = False
            if fc * f2 <= 0:
                v_g1 = c
                f1 = fc
            elif fc * f1 <= 0:
                v_g2 = c
                f2 = fc
            fvg = fabs(fc / v_a)

        return c

    def getTsat(self, p: float):
        H2O = pm.get('mp.H2O')
        Tsat = H2O.T(p=p * self.xni[1], x=0.00001)
        return Tsat  # K

    @staticmethod
    def getpsat(T: float):
        H2O = pm.get('mp.H2O')
        psat = H2O.p(T=T, x=0.00001)
        return psat  # Pa

    def Wilke_visc_calculator(self, T: float):

        def Sutherlands_law(comp: str):
            mu_0, T_0, S_mu = 0.0, 0.0, 0.0
            if comp == 'CO2':
                mu_0, T_0, S_mu = 1.370E-5, 273, 222
            elif comp == 'N2':
                mu_0, T_0, S_mu = 1.663E-5, 273, 107
            elif comp == 'O2':
                mu_0, T_0, S_mu = 1.919E-5, 273, 139
            elif comp == 'steam':
                mu_0, T_0, S_mu = 1.12E-5, 350, 1064
            mu = mu_0 * ((T / T_0) ** (3 / 2)) * (T_0 + S_mu) / (T + S_mu)
            return mu

        def Wilke_mix_rule():
            eta_frac = 0.0
            for i in range(len(self.xni)):
                den = 0.0
                for j in range(len(self.xni)):
                    mw_rel = self.mw_elem[i] / self.mw_elem[j]
                    op = self.xni[j] * (1 / sqrt(8)) * ((1 + mw_rel) ** (-0.5))
                    den += op * ((1 + (((self.mu_elem[i] / self.mu_elem[j]) ** 0.5) * (mw_rel ** 0.25))) ** 2)
                eta_frac += self.xni[i] * self.mu_elem[i] / den
            return eta_frac

        mu_CO2, mu_N2, mu_O2 = Sutherlands_law('CO2'), Sutherlands_law('N2'), Sutherlands_law('O2')
        mu_steam = Sutherlands_law('steam')
        self.mu_elem = [mu_CO2, mu_steam, mu_N2, mu_O2]
        self.mu_mix = Wilke_mix_rule()

        return self.mu_mix


def main():
    mix = mixpm("ig")
    mix.rel_error = 1E-9
    mix.setmix(10, 22, 2)
    pressure = mix.get_props_by_Tpd({'T': 1800, 'd': 2}, 'p')
    print(f'Resulta una presión de {(pressure / 101300):.7f} atm')
    enthalpy = mix.get_props_by_Tpd({'T': 1800, 'd': 2}, 'h')
    print(f'Resulta una entalpía de {enthalpy:.7f} kJ/kg')
    entropy = mix.get_props_by_Tpd({'T': 1800, 'd': 2}, 's')
    print(f'Resulta una entropía de {entropy:.7f} kJ/kgK')
    pr = mix.get_props_with_hs({'s': entropy, 'T': 1800}, {'p': 10 * 101300})
    print(f'Resulta una presión de {(pr / 101300):.7f} atm')
    Temperature = mix.get_props_with_hs({'h': enthalpy}, {'T': 1500})
    print(f'Resulta una temperatura de {Temperature:.7f} K')
    entropy2 = mix.get_props_by_Tpd({'T': Temperature, 'p': pressure}, 's')
    print(f'Resulta una entropía de {entropy2:.7f} kJ/kgK')
    return


if __name__ == '__main__':
    main()
