"""
Este módulo se encarga de simular el comportamiento de un difusor
"""

import math
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

units = {'T': 'K', 'P': 'Pa', 'H': 'J/kg'}


math.pi
pi = math.pi

def secant_method(func, x0, x1, P_out, tol=1e-6, max_iter=1000):
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
        x2 = x1 - func(x1,P_out)[0] * ((x1 - x0) / (func(x1,P_out)[0] - func(x0,P_out)[0]))
        dif = func(x2,P_out)[0]
        if dif < tol:
            return func(x2,P_out)
        x0, x1 = x1, x2
        iter_count += 1
    print("El metodo de la secante no convergio despues de", max_iter, "iteraciones.")
    return None  

def calculateMassFlow(A,v,rho): #Funcion que permite el calculo del flujo masico
            m_dot = rho * v * A
            return m_dot

def difusor(variables): #Funcion principal del modulo
    '''
    Variables de entrada:
        - H --> Altura de vuelo [m]
        - v --> Velocidad de vuelo [m/s]
        - R_e --> Radio de entrada [m]
        - rm_comp --> Radio medio del compresor [m]
        - H_comp --> Altura del primer alabe del compresor [m]
        - r_p --> relacion de presiones
    '''
    try:
        variables_difusor = {**variables[0], **variables[1]}
        H = float(variables_difusor['H_vuelo'])
        v = float(variables_difusor['v_vuelo'])
        R_e = float(variables_difusor['R_e'])
        rm_comp = float(variables_difusor['rm'])
        H_comp = float(variables_difusor['H'])
        r_p = float(variables_difusor['r_p'])
        
        #Calcular Temperatura y Presion de entrada [ecuaciones en K y Pa]
        if H > 25000:
            
            T_in = -131.21 + 0.00299*H + 273
            P_in = (2.488 * (T_in/ 216.6)**-11.388)*1000
        
        elif 11000 < H < 25000 :
            
            T_in = -56.46 + 273
            P_in = (22.65 * math.exp(1.73 - .000157 * H))*1000
        
        else:
            
            T_in = 15.04 - .00649 * H + 273
            P_in = (101.29 * ((T_in/288.08)**5.256))*1000
        
        #Entalpias J/kg
        h_a = CP.PropsSI('H', 'T', T_in + 273, 'P', P_in, 'air')
        hT_a = h_a + (v**2)/2
        hT_in = hT_a
        hT_out = hT_in 
        
        #Numero de Mach
        a =  PropsSI('A', 'T', T_in, 'P', P_in, 'air')  # Velocidad del sonido [m/s]
        M = v/a
        
        #Temperaturas
        cp = PropsSI('CPMASS', 'T', T_in, 'P', P_in, 'air')  # [J/kg*K]
        cv = PropsSI('CVMASS', 'T', T_in, 'P', P_in, 'air')  # [J/kg*K]
        gamma = cp/cv
        TT_in = T_in*(1+((gamma-1)/2)*(M**2))
        
        #Presiones [Pa]
        rho = CP.PropsSI('D', 'T', T_in, 'P', P_in, 'air') #[kg/m^3]
        PT_in = P_in*((1+((gamma-1)/2)*(M**2))**(gamma/(gamma-1)))
        
        #Areas [m^2]
        A_in = math.pi*(R_e**2)
        A_out = 2*math.pi*rm_comp*H_comp
        
        #Flujo masico [kg/s]
        m_dot = calculateMassFlow(A_in, v, rho)
        
        #propiedades totales salida
        TT_out = TT_in
        PT_out = r_p*PT_in
    
            
        def T_salida(T, P):
            rho_out = CP.PropsSI('D', 'T', T, 'P', P, 'air') #densidad de salida [kg/m^3]
            a_out = PropsSI('A', 'T', T, 'P', P, 'air') #velocidad del sonido [m/s]
            c_out = m_dot/(A_out*rho_out) #valocidad de saida [m/s]
            M_out = c_out/a_out #Mach de salida
            cp_out = PropsSI('CPMASS', 'T', T, 'P', P, 'air') #cp a la salida [J/kg*K]
            cv_out = PropsSI('CVMASS', 'T', T, 'P', P, 'air')  #cv a la salida [J/kg*K]
            gamma_out = cp_out/cv_out
            T_out = TT_out/(1+((gamma_out-1)/2)*(M_out**2)) #Temperatura a la salida [K]
            P_out = PT_out/((1+((gamma_out-1)/2)*(M_out**2))**(gamma_out/(gamma_out-1))) #Presion a la salida [Pa]
            dif = T_out-T
            return dif, T_out, c_out, P_out 
        
        dif, T_out, c_out, P_out = secant_method(T_salida, 0.9*T_in, 2.1*T_in, 0.9*P_in) # se usa el metodo de la secante para encontrar la Temperatura y Presion a la salida
        
        '''
        resultados de interes del difusor:
            T_e : Temperatura de entrada calculada para el difusor [K]
            P_e : Presion de entrada calculada para el difusor [Pa]
            Ca_1 : Velocidad de salida del difusor/ de entrada al compresor [m/s]
            T_e_comp : Temperatura de salida del difusor / de entrada al compresor [K]
            P_e_cop : Presion de salida del difusor / de entrada al compresor [Pa]
            m_dot : Fujo masico calculado para el motor
            v : velocidad de vuelo
        '''
        resultados_difusor = {
            'T_e' : T_in,
            'P_e' : P_in,
            'Ca_1' : c_out,
            'T_e_comp': T_out,
            'P_e_comp': P_out,
            'm_dot': m_dot,
            }
        
        print('resultados difusor:', resultados_difusor)
        
        
        
    except ValueError:
        raise("Error", "Ingresa valores validos en todos los campos")
    
    return resultados_difusor

if __name__ == '__main__':
    variables = [{
        'H_vuelo' : 8500,
        'v_vuelo' : 220, 
        'R_e' : 0.106,
        'r_p': 0.98,
        'k': 'None'
        }, {
        'H' : 0.04287,
        'rm': 0.21425}]
    
    i = 0
    m_dot = 0
    
    print(difusor(variables, i, m_dot))