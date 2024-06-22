import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import font
from tkinter import Frame
from tkinter import Canvas
from scipy import linalg
import math
from numpy import matrix
import numpy as np
import CoolProp.CoolProp as CP
import pandas as pd
from openpyxl import load_workbook
import openpyxl
import re
from tkinter import scrolledtext
###############################################
math.pi
pi = math.pi

def abrir_txt():
    # Abre el archivo de configuración en modo lectura
    # Abre el archivo de configuración en modo lectura
    with open('datos.txt', 'r') as archivo:
        # Lee todas las líneas del archivo
        lineas = archivo.readlines()
    
    # Inicializa un diccionario para almacenar las variables y sus valores
    variables = {}
    
    # Recorre cada línea del archivo
    for linea in lineas: 
        if not linea.startswith('#'):
            # Divide la línea en nombre de variable y valor utilizando el signo "=" como separador
            nombre, valor = linea.strip().split('=')
            # Almacena el nombre de la variable y su valor en el diccionario
            variables[nombre.strip()] = valor.strip()
        else:
            continue
    return variables

def biseccion(func, a, b, tol=1e-6, max_iter=1):
    if func(a)[0] * func(b)[0] >= 0:
        raise ValueError("La función debe cambiar de signo en el intervalo [a, b].")

    for t in range(max_iter):
        c = (a + b) / 2
        diferencia_c, T_3_c, C_3_for_c, h_3s_c, rho_3_c = func(c)
        if abs(diferencia_c) < tol:
            return c, T_3_c, C_3_for_c, h_3s_c, rho_3_c
        elif func(a)[0] * diferencia_c < 0:
            b = c
        else:
            a = c
    raise ValueError("El método de bisección no convergió después de {} iteraciones.".format(max_iter))    

def validar_nombre(nombre):
    # Utilizar una expresión regular para verificar que el nombre solo contenga letras y números
    #return re.match("^[a-zA-Z0-9]+$", nombre)
    return True 


def guardar_resultados_en_txt(resultados,j):
    try:
        # Abre el archivo en modo escritura
        with open('datos.txt', 'a') as archivo:
            # Agrega un comentario al principio del archivo
            archivo.write(f"#RESULTADOS COMPRESOR ETAPA {j}\n")

            # Itera sobre las claves y valores del diccionario
            for clave, valor in resultados.items():
                # Escribe el nombre de la variable y su valor en el archivo
                archivo.write(f'{clave} = {valor}\n')
            
            print('Etapa', j, 'guardada con éxito')
        return True
    except Exception as e:
        print(f"Error al guardar en .txt: {e}")
        return False

def calculateMasFlow(P_e,R,T_e,rm,H,Ca_1):
            Area = 2 * pi * rm * H
            rho = CP.PropsSI('D', 'T', T_e, 'P', P_e, 'air')
            m_dot = rho * Ca_1 * Area
            return m_dot
def losses(angulo_1,angulo_2,alpha_m,Altura,k,i,i_ref,S,l):
            sigma = 1/(S/l)
            if i < 0:
                messagebox.showerror("Error", "No se puede calcular Lieblein con un valor de 'i' negativo.")
            else:
                Lieblein = (math.cos(angulo_2))/(math.cos(angulo_1))*(1.12 + k*((i*180/pi) - i_ref)**1.43 + (0.61/sigma)*math.cos(angulo_1)**2* (math.tan(angulo_1)-math.tan(angulo_2)))
                if Lieblein < 0:
                    mensaje_error = "No se puede continuar: Lieblein es {Lieblein}"
                    raise ValueError(mensaje_error)
                elif Lieblein < 0:
                    return -1
                if Lieblein > 2.35:
                    print("Se ha superado el valor límite de Lieblein. El programa se detendrá.")
                    raise Exception("Valor de Lieblein superado")
            Xi_p = 2 * ((0.004)/(1-1.17*math.log(Lieblein)))*sigma*(math.cos(angulo_1)**2/math.cos(angulo_2)**3)
            Xi_a = 0.02*S/Altura *sigma *(math.cos(angulo_1)**2)/(math.cos(alpha_m)**3)
            C_L_r = 2 * 1/sigma * (math.tan(angulo_1)-math.tan(angulo_2))* math.cos(alpha_m)
            Xi_s = 0.018 * C_L_r**2 * sigma * (math.cos(angulo_1)**2)/(math.cos(alpha_m)**3)
            Xi_global = Xi_p + Xi_a + Xi_s
            #print(Xi_global)
            return Xi_global, Lieblein
def desviation(desv,beta_1_prima, beta_2_prima, S, l):
            theta = beta_1_prima - beta_2_prima #curvatura perfil
            m = 0.23 * (2*0.5)**2 + desv/500 #0.5 debido a que es un arco de circunferencia
            n = 0.5
            sigma = 1/(S/l) #solidez
            delta = m * theta * (1/sigma)**n
            return delta






def calcular(variables):
    res_c = []
    res_v = []
    j=[]
    try:
        num_etapas = int(variables['num_etapas'])
        alpha_1_grad = float(variables['alpha_1_grad'])
        beta_1_prima_gr = float(variables['beta_1_prima_gr'])
        beta_2_prima_gr = float(variables['beta_2_prima_gr'])
        alpha_2_prima_gr = float(variables['alpha_2_prima_gr'])
        alpha_3_prima_gr = float(variables['alpha_3_prima_gr'])
        rm = float(variables['rm'])
        H = float(variables['H'])
        H_2 = float(variables['H_2'])
        H_3 = float(variables['H_3'])
        S = float(variables['S'])
        l = float(variables['l'])
        Ca_1 = float(variables['Ca_1'])
        P_e = float(variables['P_e'])
        T_e = float(variables['T_e'])
        w_rpm = float(variables['w_rpm'])
        T_out_cc = float(variables['T_out'])
        T_ref = float(variables['T_ref'])
        PC = float(variables['PC'])
        eta = float(variables['eta'])
        x = float(variables['x'])
        y = float(variables['y'])
        r_p = float(variables['r_p'])
        cte_H_2 = H_2/H
        cte_H_3 = H_3/H_2 
        Rm = float(variables['Rm'])
        H_1_t = float(variables['H_1'])
       

        
        
        for etapa in range(1, num_etapas + 1):
            try:
                print("Iniciando cálculos para etapa", etapa)
                
                if etapa > 1:
                    print("Actualizando valores para etapa", etapa)
                    P_e = P_3
                    T_e = T_3
                    Ca_1=Ca_3
                    alpha_1 = alpha_3
                    H = H_3
                    H_2 = cte_H_2 * H
                    H_3 = cte_H_3 * H_2
                    print("Valores asignados:", P_e, T_e, Ca_1, alpha_1, H, H_2, H_3)   
                    
                        
                ###############################
                if etapa == 1:
                     alpha_1 = alpha_1_grad/180 * pi #Esto se hace para que no cambie luego después de las siguientes etapas

                #alpha_1 = alpha_1_grad/180 * pi #Angulo real del flujo
                beta_1_prima = beta_1_prima_gr/180 * pi #Angulo del álabe a la entrada (ETAPA POR ETAPA)
                beta_2_prima = beta_2_prima_gr/180  * pi #Angulo salida del rotor relativo (ángulo álabe salida rotor) (ETAPA POR ETAPA)
                alpha_2_prima = alpha_2_prima_gr/180 * pi #Angulo del alabe del estator (falta por poner)
                alpha_3_prima = alpha_3_prima_gr/180 * pi

                k = 0.0117 #cte para álabes NACA-65
                i_ref = 0/180 * pi
                R = 287
                rho = CP.PropsSI('D', 'T', T_e, 'P', P_e, 'air')
                mas_flow = calculateMasFlow (P_e,R,T_e,rm,H,Ca_1)
                print('El flujo másico es:', mas_flow)
                #CÁLCULO DE VELOCIDAD AXIAL MÁXIMA
                beta_1_max = beta_1_prima
                U = w_rpm *(pi/30) * rm
                #print(U)

                # Imprimir la solución
                #CALCULO DE DESVIACION
                beta_1_prima = beta_1_max
                alpha_rotor = beta_2_prima*180/pi
                delta_rotor = desviation(alpha_rotor,beta_1_prima, beta_2_prima, S, l)
                delta_rotor_rad = delta_rotor*pi/180
                beta_2 = beta_2_prima + delta_rotor_rad
                
                #TRIANGULO DE VELOCIDADES
                C_1 = Ca_1/ math.cos(alpha_1)
                U = w_rpm *(pi/30) * rm
                C_u1= C_1 * math.sin(alpha_1)
                W_u1= U - C_u1
                W_1 = math.sqrt(W_u1**2 + Ca_1**2)
                beta_1_rad = math.atan(W_u1/Ca_1) #Angulo relativo del flujo a la entrada del álabe
                i=beta_1_rad - beta_1_prima
                i_grad = i*180/pi
                beta_1 = beta_1_rad*180/pi
                theta_1_rad = 180/180 * pi - (90/180 * pi + beta_2)
                theta_1 = theta_1_rad*180/pi
                h_1 = CP.PropsSI('H', 'T', T_e, 'P', P_e, 'air')
                h_01 = h_1 + (C_1**2)/2 #J/Kg
                a_1 = CP.PropsSI('A', 'T', T_e, 'P', P_e, 'air')
                M_1 = C_1/a_1
                M_1_relativo = W_1/a_1
                cp_1 = CP.PropsSI('C', 'T', T_e, 'P', P_e, 'air') 
                print(cp_1)
                cv_1 = CP.PropsSI('CVMASS', 'T', T_e, 'P', P_e, 'air') 
                gamma_1 = cp_1/cv_1

                T_01d = T_e * (1 + (gamma_1 - 1)/2 * M_1**2)#Comprobación mediante fórmula (no se va a usar)
                rho_1 = CP.PropsSI('D', 'T', T_e, 'P', P_e, 'air')

                s_1 = CP.PropsSI('S', 'P', P_e, 'T', T_e, 'air')
                T_01 = CP.PropsSI('T', 'H', h_01, 'S', s_1, 'air') #Esta es mejor ya que usamos las tablas
                P_01 = CP.PropsSI('P', 'T', T_01, 'S', s_1, 'air')
                s_01 = CP.PropsSI('S', 'P', P_01, 'T', T_01, 'air') #Tienen que dar lo mismo que s_1
                H_1 = mas_flow/(rho_1 * Ca_1 * 2 * pi * rm)                  

                def funcion(Ca_2):
                            # Función diferencia_2 = (m_dot_for_2 - mas_flow)
                            # El resultado de esta función es lo que queremos que cambie de signo (se haga cero).
                            # Aquí debes incluir el cálculo de la diferencia y devolver el resultado.
                            # Por ejemplo:
                            W_2 = Ca_2 * (math.sin(90/180 * pi)/ math.sin(theta_1_rad))
                            W_u2 = W_2 * math.sin(beta_2)
                            X = (math.sin(beta_1_rad)*W_1)/(math.sin(90/180 * pi))- (W_u2)
                            C_u2 = X + C_u1
                            C_2 = math.sqrt(Ca_2**2 + C_u2**2)
                            alpha_2 = math.atan(C_u2/Ca_2) #Angulo real salida rotor y coincide con la entrada al estator
                            alpha_2_grad = alpha_2 * 180/pi
                            W = U * (C_u2-C_u1)/1000
                            Cp_1 = CP.PropsSI('C', 'T', T_e, 'P', P_e, 'air') 
                            h_02 = W*1000 + h_01
                            h_02 = float(h_02)
                            T_02 = ((h_02 - h_01)/Cp_1)+ T_01 
                            T_02 =float(T_02)
                            s_02s = s_01
                            s_2s = s_1
                            s_2s =float(s_2s)
                            #Calculamos las perdidas del rotor
                            angulo_1 = beta_1_rad 
                            angulo_2 = beta_2
                            alpha_m_rotor = (beta_1_rad + beta_2)/2

                            Altura = H_2 #Valor conocido
                            Xi_rotor, Lieblein_rotor = losses(angulo_1,angulo_2,alpha_m_rotor,Altura,k,i,i_ref,S,l)

                            h_2 = h_02 - ((C_2**2)/2)
                            T_2 = T_02 - (C_2** 2)/(2*Cp_1)#COOLPROP
                            h_2s = h_2 - Xi_rotor * (C_2**2/2)
                            h_2s =float(h_2s)
                            P_2 = CP.PropsSI('P', 'H', h_2s, 'S', s_2s, 'air')
                            P_2 =float(P_2)
                            rho_2 = CP.PropsSI('D', 'T', T_2, 'P', P_2, 'air')
                            m_dot_for_2 = rho_2 * C_2 * math.cos(alpha_2) * 2 * np.pi * rm * H_2
                            diferencia_2 = m_dot_for_2 - mas_flow
                            return diferencia_2, h_2, C_u2, C_2, T_02

                # Estimaciones iniciales del intervalo [a, b] donde se encuentra la raíz.
                a = 1
                b = 300
                print(a)
                # Llamada al método de bisección para encontrar la raíz
                Ca_2, h_2, C_u2, C_2, T_02 = biseccion(funcion, a, b)
                print(b)

                W = U * (C_u2-C_u1)/1000 #Las unidades son KJ/Kg
                print ('W: ', W)
                W_dot = W * mas_flow 
                print (W_dot) #KW

                Cp_1 = CP.PropsSI('C', 'T', T_e, 'P', P_e, 'air')
                #Calculamos las perdidas del rotor
                angulo_1 = beta_1_rad 
                angulo_2 = beta_2
                Altura = H_2 #Valor conocido
                alpha_m_rotor = (beta_1_rad + beta_2)/2
                alpha_m_rotor_grad = alpha_m_rotor*180/pi
                Xi_rotor, Lieblein_rotor = losses(angulo_1,angulo_2,alpha_m_rotor,Altura,k,i,i_ref,S,l)
                s_02s = s_01
                s_2s = s_1
                alpha_2 = math.atan(C_u2/Ca_2)

########################################################################
                T_2 = T_02 - (C_2 ** 2)/(2*Cp_1)
                h_02 = h_2 + (C_2**2/2)
                h_02 =float(h_02)
                h_2s_nuevo = h_2 - Xi_rotor * (W_1**2/2)
                P_2s = CP.PropsSI('P', 'H', h_2s_nuevo, 'S', s_2s, 'air')
                P_2 = P_2s
                P_2 = float(P_2)
                rho_2 = mas_flow/(H_2 * C_2 * math.cos(alpha_2)* 2 * pi * rm)
                s_2 = CP.PropsSI('S', 'H', h_2, 'P', P_2, 'air')
                s_2 = float(s_2)
                s_02 = s_2
                s_02 = float(s_02)
                P_02 = CP.PropsSI('P', 'H', h_02, 'S', s_02, 'air')
                h_02s = CP.PropsSI('H', 'P', P_02, 'S', s_1, 'air')
                T_2s = CP.PropsSI('T', 'P', P_2, 'S', s_2, 'air')#Esto es para comprobar T_2s = T_2
                print('T_2',T_2s)
                #CALCULO DE DESVIACION ALPHA_3
                alpha_estator = alpha_3_prima*180/pi
                delta_estator = desviation(alpha_estator,beta_1_prima, beta_2_prima, S, l)
                delta_estator_rad = delta_estator*pi/180
                alpha_3 = alpha_3_prima + delta_estator_rad
                h_03 = h_02
                #Calculamos las perdidas del estator
                angulo_1 = alpha_2
                angulo_2 = alpha_3
                Altura_estator = H_3
                print(Altura_estator)
                i_estator = alpha_2 - alpha_2_prima
                i_estator_grad = i_estator*180/pi
                alpha_m_estator =(alpha_2 + alpha_3)/2
                alpha_m_estator_grad = alpha_m_estator*180/pi
                Xi_estator, Lieblein_estator = losses(angulo_1,angulo_2,alpha_m_estator,Altura_estator,k,i_estator,i_ref,S,l)
                
                Cp_2 = CP.PropsSI('C', 'T', T_2, 'P', P_2, 'air')
                cp_1
                print(Cp_2)
                T_03_for = ((h_03 - h_02)/Cp_2)+ T_02 
                print('T3',T_03_for)
                def funtion(h_3_ini):

                            C_3_for = math.sqrt(2*(h_03 - h_3_ini))
                            T_3 = T_03_for - (C_3_for** 2)/(2*Cp_2)
                            h_3s = h_3_ini - Xi_estator * (C_3_for**2/2)
                            P_3 = CP.PropsSI('P', 'H', h_3s, 'S', s_2, 'air')
                            P_3 =float(P_3)
                            s_3 = CP.PropsSI('S', 'P', P_3, 'T', T_3, 'air')
                            P_03 = CP.PropsSI('P', 'T', T_03_for, 'S', s_3, 'air') 
                            rho_3 = CP.PropsSI('D', 'T', T_3, 'P', P_3, 'air')
                            C_3_axial = C_3_for * math.cos(alpha_3)
                            m_dot_for = rho_3 * C_3_axial *2 * pi * rm * H_3
                            diferencia = (m_dot_for - mas_flow)
                            return diferencia, T_3, C_3_for, h_3s, rho_3

                
                # Estimaciones iniciales del intervalo [a, b] donde se encuentra la raíz.
                a = h_2
                b = h_03

                # Llamada al método de bisección para encontrar la raíz
                h_3_ini, T_3, C_3_for, h_3s, rho_3 = biseccion(funtion, a, b)

                C_3 = math.sqrt(2*(h_03 - h_3_ini))
                T_3 = T_03_for - (C_3** 2)/(2*Cp_2)
                h_3s = h_3_ini - Xi_estator * (C_3**2/2)
                h_3s =float(h_3s)
                s_2 = float(s_2)
                P_3 = CP.PropsSI('P', 'H', h_3s, 'S', s_2, 'air')
                P_3 =float(P_3)
                rho_3 = CP.PropsSI('D', 'T', T_3, 'P', P_3, 'air')
                s_3 = CP.PropsSI('S', 'P', P_3, 'T', T_3, 'air') #Da igual poner P_3 que P_3_nuevo, dan lo mismo
                P_03 = CP.PropsSI('P', 'T', T_03_for, 'S', s_3, 'air') 
                T_03ss = CP.PropsSI('T', 'P', P_03, 'S', s_1, 'air')
                T_3ss = CP.PropsSI('T', 'P', P_3, 'S', s_1, 'air')
                h_3ss = CP.PropsSI('H', 'T', T_3ss, 'S', s_1, 'air')
                h_03ss = CP.PropsSI('H', 'T', T_03ss, 'S', s_1, 'air')#Da igual poner entropia, que temperatura o presion, va a dar lo mismo
                h_03s = h_02
                rp = P_3/P_e
                Ca_3 = C_3 * math.cos(alpha_3) 
                #RENDIMIENTOS
                eta_TT_2 = (h_03ss - h_01)/(h_03 - h_01)
                eta_TT_3 = 1 - ((h_3_ini - h_3s) + (T_3/T_2)*(h_2 - h_2s_nuevo))/(h_03 - h_01) 
                eta_TT_4 = 1 - ((h_3_ini - h_3ss)/(h_03 - h_01))
                               
                
                resultados={
                'Etapa': etapa,
                'U': w_rpm,
                'H':H,
                'H_2':H_2,
                'H_3':H_3,
                'Ca_1': Ca_1,
                'C_1': C_1,
                'P_e': P_e,
                'T_e': T_e,
                'h_1': h_1,
                'h_01': h_01,
                's_1': s_1,
                'rho_1': rho_1,
                'Cp_1': Cp_1,
                'Ca_2': Ca_2,
                'C_2': C_2,
                'P_2': P_2,
                'T_2': T_2,
                'h_2': h_2,
                'h_02': h_02,
                's_2': s_2,
                'rho_2': rho_2,
                'Cp_2': Cp_2,
                'Ca_3': Ca_3,
                'C_3': C_3,
                'P_3': P_3,
                'T_3': T_3,
                'h_3_ini': h_3_ini,
                'h_03': h_03,
                's_3': s_3,
                'rho_3': rho_3,
                'm_dot': mas_flow,
                'eta_TT': eta_TT_4,
                'Rp': rp,
                'Pérdida rotor': Xi_rotor,
                'Pérdida estátor': Xi_estator,
                'Potencia': W_dot,
                'M_1 relativo': M_1_relativo,
                'Lieblein rotor': Lieblein_rotor,
                'Lieblein estator': Lieblein_estator,
                'alpha_1': alpha_1,
                'beta_1_rad': beta_1_rad,
                'beta_1_prima': beta_1_prima,
                'beta_2': beta_2,
                'beta_2_prima': beta_2_prima,
                'alpha_2': alpha_2,
                'alpha_2_prima': alpha_2_prima,
                'alpha_3': alpha_3,
                'alpha_3_prima': alpha_3_prima,
                'P_03' : P_03,
                'T_out_cc' : T_out_cc,
                'T_ref' : T_ref,
                'PC' : PC,
                'eta' : eta,
                'x' : x,
                'y' : y,
                'r_p' : r_p,
                'Rm': Rm,
                'H_1_t': H_1_t
                }
                print(resultados)
                resultados_clave = list(resultados)
                resultados_valor = list(resultados.values())
                j.append(etapa)
                res_c.append(resultados_clave)
                res_v.append(resultados_valor)
                
               
                

            except ValueError:
                print("Error en cálculos de etapa", etapa)
                messagebox.showerror("Error", "Ingresa valores numéricos válidos en los cálculos de esta etapa")
            
            
        print("Cálculos realizados correctamente")
        
    except ValueError:
        
        print("Error", "Ingresa valores válidos en todos los campos")
    
    print('listo')
    return resultados   
        


if __name__ == '__main__':
    calcular(abrir_txt())

