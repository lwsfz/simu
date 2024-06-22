"""
Modulo encargado de simular el comportamiento de la tobera
"""

import math
import CoolProp
from CoolProp.CoolProp import PropsSI

def main_tobera(variables, resultados_dif, resultados_turb, composicion_cc):
    try:
        
        '''
        eta_n --> Rendimiento de la tobera
        tobera --> Geometria de la tobera
        R_in --> Radio de entrada de la tobera [m]
        H_in --> Altura a la salida del ultimo escalonamiento de la turbina [m]
        P_ext --> Presion exterior [Pa]
        PT_in --> Presion total de entrada a la tobera [Pa]
        TT_in --> Temperatura total de entrada a la tobera [K]
        m_dot --> Flujo masico [kg/s]
        n_CO2 --> kmolCO2/kmol_cble
        n_H2O --> kmolH2O/kmol_cble
        n_N2 --> kmolN2/kmol_cble
        n_O2_h --> kmolO2(humos)/kmol_cble
        '''
        
        variables_tobera =  {**variables[4], **resultados_dif, **resultados_turb, **composicion_cc}
        eta_n = float(variables_tobera['eta_n'])
        tobera = str(variables_tobera['tobera'])
        R_in = float(variables_tobera['R_in'])
        H_in = float(variables_tobera['H_in'])
        P_ext = float(variables_tobera['P_e'])
        PT_in = float(variables_tobera['PT_out'])
        TT_in = float(variables_tobera['TT_out'])
        m_dot = float(variables_tobera['m_dot'])
        n_CO2= float(variables_tobera['n_CO2'])
        n_H2O = float(variables_tobera['n_H2O'])
        n_N2 = float(variables_tobera['n_N2'])
        n_O2_h = float(variables_tobera['eta_n'])
        M_cble = float(variables_tobera['M_cble'])
        
        #Masas molares kg/kmol
        M_CO2 = 44 
        M_N2 = 28 
        M_H2O = 18
        M_O2 = 32 
        
        #Funcion para calular la relacion de masa_i/masa_cble
        def frac(n,M): #n = moles_i/moles_cble y M = masa molar_i
            m = n*M/M_cble #kg_i/kg_cble
            return m  
        
        #funcion para obetener la fraccion masica de los productos
        def fm_productos(n,M):
            m_flujo = frac(n_CO2,M_CO2)+frac(n_N2,M_N2) + frac(n_O2_h,M_O2) + frac(n_H2O,M_H2O)
            fm = frac(n,M)/m_flujo
            return fm
       
        cp_in = fm_productos(n_CO2,M_CO2)*PropsSI('CPMASS', 'T', TT_in, 'P', PT_in, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('CPMASS', 'T', TT_in, 'P', PT_in, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('CPMASS', 'T', TT_in, 'P', PT_in, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('CPMASS', 'T', TT_in, 'P', PT_in, 'Nitrogen')  # [J/kg*K]
        cv_in = fm_productos(n_CO2,M_CO2)*PropsSI('CVMASS', 'T', TT_in, 'P', PT_in, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('CVMASS', 'T', TT_in, 'P', PT_in, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('CVMASS', 'T', TT_in, 'P', PT_in, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('CVMASS', 'T', TT_in, 'P', PT_in, 'Nitrogen')  # [J/kg*K]
        gamma = cp_in/cv_in
        
        
        '''
        Se definen a continuacion los funcionamiento de las toberas
        '''
        if tobera == 'convergente':
            M_out = 1
            #Calculo P_out
            P_out = PT_in*(1/eta_n*(1/(1+(gamma-1)/2*M_out)-1+eta_n))**(gamma/(gamma-1))
            
            if P_out < P_ext: #Si esto ocurre significa que M<1 y la suposicion inicial esta mal
                P_out = P_ext #Se supone tobera adaptada
                
                #se calcula Mach
                M_out = math.sqrt(2/(gamma-1)*(1/((P_out/PT_in)((gamma-1)/gamma)*eta_n+1-eta_n)-1))
                
            else: 
                pass
                
            
            #se calcula T_out,s (isoentropica) para posteriormente calcular T_out
            T_out_s = TT_in*((P_out/PT_in)**((gamma-1)/gamma))
            h_out_s = fm_productos(n_CO2,M_CO2)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            h_Tin = fm_productos(n_CO2,M_CO2)*PropsSI('H', 'T', TT_in, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('H', 'T', TT_in, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('H', 'T', TT_in, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('H', 'T', TT_in, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            h_out = h_Tin - eta_n*(h_Tin-h_out_s)
            c_out = math.sqrt(2*(h_Tin-h_out))
            T_out = fm_productos(n_CO2,M_CO2)*PropsSI('T', 'H', h_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('T', 'H', h_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('T', 'H', h_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('T', 'H', h_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            
            #Calcular cp_out
            cp_out = fm_productos(n_CO2,M_CO2)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            cv_out = fm_productos(n_CO2,M_CO2)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            gamma_out = cp_out/cv_out
            rho_out = fm_productos(n_CO2,M_CO2)*PropsSI('D', 'T', T_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('D', 'T', T_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('D', 'T', T_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('D', 'T', T_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            Area_out = m_dot/(rho_out*c_out)
            Area_in = 2*math.pi*R_in*H_in
            
            if Area_out > Area_in or Area_out == Area_in:
                raise('El area de la garganta de la tobera no puede ser mayor que el area de entrada de la misma')
            
            else:
                pass
            
            '''
            Resultados de interes de la tobera:
                tobera: tipo de tobera que se esta usando (C o C-D)
                m_dot: FLujo masico del motor
                M_out: Mach a la salida de la tobera
                c_out: Velocidad a la salida de la tobera
                P_out: Presion a la salida de la tobera
                P_e: Presion en el exterior
                A_out: Area salida de la tobera
            '''
            
            resultados_tob = {
                'tobera': tobera,
                'm_dot': m_dot,
                'M_out': M_out,
                'c_out': c_out,
                'P_out': P_out,
                'P_e': P_ext,
                'A_out': Area_out
                }
        
        elif tobera == 'convergente-divergente':
            
            #Condicion de tobera adaptada para calcular el Mach a la salida
            P_out = P_ext
            M_out = math.sqrt(2/(gamma-1)*(1/((P_out/PT_in)**((gamma-1)/gamma)*eta_n+1-eta_n)-1))
            
            #Con el Mach obtenido se determina la relacion de areas A_out/A*
            def get_ra(M):
                A = (1/M)*((2/(gamma+1))*(1+(gamma-1)*(M**2)/2))**(-1/2)
                B = eta_n*(1+(1-gamma)/(eta_n*(1+gamma)))
                C = 1/(1+(gamma-1)*M**2/2) - 1 +eta_n
                rel_area = A*((B/C)**(gamma/(gamma-1)))
                return rel_area
            
        
            rel_area = get_ra(M_out)
            
            #se calcula T_out,s (isoentropica) para posteriormente calcular T_out
            #se calcula T_out,s (isoentropica) para posteriormente calcular T_out
            T_out_s = TT_in*((P_out/PT_in)**((gamma-1)/gamma))
            h_out_s = fm_productos(n_CO2,M_CO2)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('H', 'T', T_out_s, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            h_Tin = fm_productos(n_CO2,M_CO2)*PropsSI('H', 'T', TT_in, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('H', 'T', TT_in, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('H', 'T', TT_in, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('H', 'T', TT_in, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            h_out = h_Tin - eta_n*(h_Tin-h_out_s)
            c_out = math.sqrt(2*(h_Tin-h_out))
            T_out = fm_productos(n_CO2,M_CO2)*PropsSI('T', 'H', h_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('T', 'H', h_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('T', 'H', h_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('T', 'H', h_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            
            #Calcular cp_out
            cp_out = fm_productos(n_CO2,M_CO2)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('CPMASS', 'T', T_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            cv_out = fm_productos(n_CO2,M_CO2)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('CVMASS', 'T', T_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            gamma_out = cp_out/cv_out
            rho_out = fm_productos(n_CO2,M_CO2)*PropsSI('D', 'T', T_out, 'P', P_out, 'CO2')+fm_productos(n_H2O,M_H2O)*PropsSI('D', 'T', T_out, 'P', P_out, 'Water')+fm_productos(n_O2_h,M_O2)*PropsSI('D', 'T', T_out, 'P', P_out, 'Oxygen')+fm_productos(n_N2,M_N2)*PropsSI('D', 'T', T_out, 'P', P_out, 'Nitrogen')  # [J/kg*K]
            Area_out = m_dot/(rho_out*c_out)
            Area_g = Area_out/rel_area
            Area_in = 2*math.pi*R_in*H_in
            
            if Area_g > Area_in or Area_g == Area_in:
                raise('El area de la garganta de la tobera no puede ser mayor que el area de entrada de la misma')
            
            else:
                pass

            #Mach despues de la onda de choque
            M_y = math.sqrt(((M_out**2)+2/(gamma_out - 1))/((2*gamma_out/(gamma_out - 1))*(M_out**2 - 1)))
        
            #Presion posterior a la onda de choque
            P_y = P_out*((1+gamma_out*M_out**2)/(1+gamma_out*M_y**2))
            print(P_y,P_ext,P_out,'presiones')
            
            if P_y < P_ext:
                raise('Existe una onda de choque en el interior de la tobera. Condicion de funcionamiento poco optima. Introduzca datos de tobera nuevamente.')
            else:
                pass
            
            '''
            Resultados de interes de la tobera:
                tobera: tipo de tobera que se esta usando (C o C-D)
                m_dot: FLujo masico del motor
                M_out: Mach a la salida de la tobera
                c_out: Velocidad a la salida de la tobera
                P_out: Presion a la salida de la tobera
                P_e: Presion en el exterior
                A_out: Area salida de la tobera
                A_g: Area de la garganta en tobera conv-div
            '''
            resultados_tob = {
                'tobera': tobera,
                'm_dot': m_dot,
                'M_out': M_out,
                'c_out': c_out,
                'P_out': P_out,
                'P_e': P_ext,
                'A_out': Area_out,
                'A_g': Area_g
                }
        
        else:
            raise(ValueError, 'Ingresa un valor str que determine la geometria de la tobera: convergente o convergente-divergente.')
            
        
    except ValueError:
        raise("Error", "Ingresa valores validos en todos los campos")
        
    return resultados_tob

if __name__ == '__main__':
    
    variables = [0,0,0,0,{
        'eta_n': 0.96,
        'tobera': 'convergente-divergente',
        'R_in': 0.1195,
        'H_in': 0.066,
        'ra': 1.05
        }]
    
    rdif = {
        'P_e': 33110.32680055705}
    
    resultados_turb = {
        'PT_out':  499261.542744159761,
        'P_out':  493103.06616336183,
        'TT_out': 1278.3459130984536,
        'c_out': 125.25050974535787,
        'm_dot': 3.8481904283202533
        }
    composicion_cc = {
        'n_CO2': 11,
        'n_H2O': 11,
        'n_N2':149.3579874,
        'n_O2_h': 39.70275614,
        'M_cble': 142
        }
    print(main_tobera(variables, rdif, resultados_turb, composicion_cc))