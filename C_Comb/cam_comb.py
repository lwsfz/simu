import pyromat as pm
import math

#Configuración de laas unidades
pm.config['unit_pressure'] = 'Pa'
pm.config['unit_temperature'] = 'K'
pm.config['unit_energy'] = 'kJ'
pm.config['unit_mass'] = 'kg'
pm.config['unit_volume'] = 'm3'

def secant_method(func, x0, x1, max_iter=10000):
    """
    Encuentra una raiz de la funcion 'func' utilizando el metodo de la secante.
    
    Args:
        func: La funcion para la cual se busca la raiz.
        x0: Primer punto inicial.
        x1: Segundo punto inicial.
        tol: Tolerancia para el criterio de parada (por defecto es 1e-6).
        max_iter: Numero maximo de iteraciones permitidas (por defecto es 1000).
        
    Returns:
        La aproximacion de la raiz encontrada.
    """
    iter_count = 0
    while iter_count < max_iter:
        x2 = x1 - func(x1)[0] * ((x1 - x0) / (func(x1)[0] - func(x0)[0]))
        dif = func(x2)[0]
        print(dif,'diferenica de exceso de aire')
        if dif < 0.00001:
            return func(x2)
        x0, x1 = x1, x2
        iter_count += 1
    print("El metodo de la secante no convergio despues de", max_iter, "iteraciones.")
    return None
    

def cam_comb(variables, resultados_comp, T_out):
    try:
        """    
        T_in,T_out,T_ref, PC, eta, x, y, P_T_in, r_p, Rm_tur, H_1
        Datos que se conocen 
        -Combustible CxHy
        -PC = Poder calorifico en KJ/kg
        -eta = rendimiento del PC
        -T_ref = Temperatura de referencia en K
        -T_in = Temperatura de entrada en K 
        -T_out = Temperatura de salida en K
        se quiere calcular el flujo de combustible a traves de la expresion del balance de energia,
        para poder conocer el valor de n,
        para ello sera necesario asumir que la combustion es completa y que 
        el calor transferido es 0.
        Una vez hechas estas suposiciones necesitamos calcular la masa del combustible y 
        los productos a traves de un balance de reacciones
        De este balance tenemos la relacion del combustible con los productos, que se 
        relacionan con las variables x e y (el H20 no es un humo no se tiene en cuenta)
        Todos los moles estan en relacion Kmol_i/Kmol_cble
        """
        variables_cam_comb = {**variables[2],**variables[3], **resultados_comp} 
         
        #Cargar datos desde el archivo txt
        T_in = float(variables_cam_comb['T_3'])
        T_ref = float(variables_cam_comb['T_ref'])
        PC = float(variables_cam_comb['PC'])
        eta = float(variables_cam_comb['eta'])
        x = float(variables_cam_comb['x'])
        y = float(variables_cam_comb['y'])
        P_T_in = float(variables_cam_comb['P_03'])
        r_p = float(variables_cam_comb['r_p'])
        m_dot = float(variables_cam_comb['m_dot'])
        Rm_tur = float(variables_cam_comb['Rm'])
        H_1 = float(variables_cam_comb['H_1'])
        
        #pyromat para propiedades
        mp_N2 = pm.get('ig.N2')
        mp_CO2 = pm.get('ig.CO2')
        mp_H2O = pm.get('ig.H2O')
        mp_O2 = pm.get('ig.O2')
        
    
        def calculate_n(n): #Funcion utilizada para el calculo del exceso de aire (n)
            #kmoles_productos/kmoles_cble
            n_CO2=x
            n_H2O = y/2 
            n_N2 = n*(x+y/4)*79/21 
            n_O2_h= (n-1)*(x+y/4)
            P_H2O = n_H2O*P_T_in/(n_CO2+n_H2O+n_N2+n_O2_h)
            
            #kmoles_reactivos/kmoles_cble
            n_O2 = n*(x + y/4)
            #n_N2 = productos
            
            #Masas molares kg/kmol
            M_CO2 = 44 
            M_N2 = 28 
            M_H2O = 18
            M_O2 = 32 
            M_cble = 12*x + y 
            
            
            #entalpias productos kJ/Kg
            h_N2_h = mp_N2.h(T_out)
            print(h_N2_h, 'h')
            h_N2_h_ref = mp_N2.h(T_ref)
            h_CO2 = mp_CO2.h(T_out)
            h_CO2_ref = mp_CO2.h(T_ref)
            
            h_H2O = mp_H2O.h(T=T_out, p=P_H2O)
            h_H2O_ref = mp_H2O.h(T=T_ref, p=101325)
            h_O2_h = mp_O2.h(T_out)
            h_O2_h_ref = mp_O2.h(T_ref)
            
            #entalpias reactivos kJ/kg
            h_N2 = mp_N2.h(T_in)
            h_N2_ref = mp_N2.h(T_ref)
            h_O2 = mp_O2.h(T_in)
            h_O2_ref = mp_O2.h(T_ref)
            
            
            '''
            CALCULO DE DOSADO
            '''
            
            #Funcion para calular la relacion de masa_i/masa_cble
            def frac(n,M): #n = moles_i/moles_cble y M = masa molar_i
                m = n*M/M_cble #kg_i/kg_cble
                return m  
            
            #Funcion para calcular la energia generada por cada producto
            def energia(n,M,h_i,h_ref):
                p = frac(n,M)*(h_i-h_ref) #kJ/kg
                return p
            
             
            #Funcion para el calculo de la fraccion masica de los reactivos
            def fm(n,M):
                m_air = frac(n_N2,M_N2) + frac(n_O2,M_O2)
                fm = frac(n,M)/m_air #kg_i/kg_air
                return fm
            
            #Entalpia del aire kJ/kg
            h_air = fm(n_N2,M_N2)*h_N2+fm(n_O2,M_O2)*h_O2 
            h_air_ref = fm(n_N2,M_N2)*h_N2_ref+fm(n_O2,M_O2)*h_O2_ref
            
            #energia aportada por el aire kJ/kg
            energia_air = h_air - h_air_ref 
            
            #energia de los productos kJ/kg
            Energia_productos = energia(n_N2,M_N2,h_N2_h,h_N2_h_ref)+energia(n_CO2,M_CO2,h_CO2,h_CO2_ref)+energia(n_O2_h,M_O2,h_O2_h,h_O2_h_ref)+energia(n_H2O,M_H2O,h_H2O,h_H2O_ref)
            
            #dosado y masa de combustible
            f = energia_air/(Energia_productos-PC*eta)
            print(f)
            m_cble = f*m_dot #kg/s
            
            #para la iteracion de n
            a = (12*x + y)/(f*(32 + 28*79/21))
            a_e = x + y/4
            n1 = a/a_e
            dif_n = n1-n
            
            "Propiedades termodinamicas necesarias para obtener c y P"
            #funcion para obetener la fraccion masica de los porductos
            def fm_productos(n,M):
                m_flujo = frac(n_CO2,M_CO2)+frac(n_N2,M_N2) + frac(n_O2_h,M_O2) + frac(n_H2O,M_H2O)
                fm = frac(n,M)/m_flujo
                return fm
            
            #Popiedades de los productos a la salida de la camcomb
            R_N2 = mp_N2.R(T_out)
            gamma_N2 = mp_N2.cp(T_out) / mp_N2.cv(T_out)
            R_O2 = mp_O2.R(T_out)
            gamma_O2 = mp_O2.cp(T_out) / mp_O2.cv(T_out)
            R_H2O = mp_H2O.R(T_out)
            gamma_H2O = mp_H2O.cp(T_out) / mp_H2O.cv(T_out)
            R_CO2 = mp_CO2.R(T_out)
            gamma_CO2 = mp_CO2.cp(T_out) / mp_CO2.cv(T_out)
            R = (fm_productos(n_N2,M_N2)*R_N2+fm_productos(n_O2_h,M_O2)*R_O2+fm_productos(n_CO2,M_CO2)*R_CO2+fm_productos(n_H2O,M_H2O)*R_H2O)*1000
            gamma = (fm_productos(n_N2,M_N2)*gamma_N2+fm_productos(n_O2_h,M_O2)*gamma_O2+fm_productos(n_CO2,M_CO2)*gamma_CO2+fm_productos(n_H2O,M_H2O)*gamma_H2O)
            
            "Proceso iterativo para obtener c y P"
            P_T_out = r_p*P_T_in
            A_in_tur = 2*math.pi*Rm_tur*H_1
            
            def get_c_out(c_out):
                P_out=P_T_out/((1+(gamma-1)/2*(c_out**2/(gamma*R*T_out)))**(gamma/(gamma-1)))
                c_out_1=m_dot*R*T_out/(P_out*A_in_tur)
                dif = abs(c_out_1 - c_out)
                c_out = c_out_1
                return dif, c_out, P_out
            
            dif, c_out, P_out = secant_method(get_c_out, 40, 200)
            
            resultados_cam_comb = { 
                'c_out' : c_out,
                'p_inlet' : P_out,
                'T_inlet' : T_out,
                'mass_flow' : m_dot,
                'n': n1,
                'm_cble': m_cble
                }
            
            composicion_cc = {
                'n_CO2': n_CO2,
                'n_H2O': n_H2O,
                'n_N2' :  n_N2,
                'n_O2_h' : n_O2,
                'M_cble': M_cble
                } 
            
            return dif_n, resultados_cam_comb, composicion_cc, n1
        
        a = 1
        b = 25
        
        dif_n, resultados_cam_comb, composicion_cc, n1 = secant_method(calculate_n,a,b)
        print(resultados_cam_comb)
        
        return resultados_cam_comb, composicion_cc
    
    except ValueError:
        print("Error", "Ingresa valores válidos en todos los campos")
    
    
    
if __name__ == '__main__':
    variables = [0,0,{
        'T_ref' : 298,
        'PC' : 40000,
        'eta' : 0.98,
        'x' : 10,
        'y' : 22
        },
        {'Rm' : 0.105,
        'H_1' : 0.055,
            }]
    
    resultados_comp = {
        'T_3' : 569.8773086415234,
        'P_03' : 642509.5219117957,
        'r_p' : 0.98,
        'm_dot' : 3.8481904283202533,
        'U' : 13000,
    
        }
    
    T_out_cc = 1392.4537895750263
    print(cam_comb(variables, resultados_comp, T_out_cc))
