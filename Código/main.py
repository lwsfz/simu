# -*- coding: utf-8 -*-
"""Funcion main del programa de simulacion que se va a usar
para llamar a todos los mÃ³dulos y codigos necesarios para 
conseguir el fin determinado"""

"Primero importamos todos los codigos"
import os
from Difusor.difusor import *
from Compresor.Compresor import * 
from C_Comb.cam_comb import *
from Turbina.solver_manager import *
from Tobera.tobera import *

"""
Funciones para abrir y guardar datos en el txt
"""

def abrir_txt(): #Abre el archivo de datos generales para el turbojet 'datos.txt'
    # Ruta a la carpeta donde está el archivo
    carpeta = os.path.join('..', 'datos')
    nombre_archivo = 'datos_supersónico.txt'

    # Ruta completa al archivo
    ruta_archivo = os.path.join(carpeta, nombre_archivo)

    with open(ruta_archivo, 'r') as archivo:
        # Lee todas las líneas del archivo
        lineas = archivo.readlines()
    
    # Inicializa diccionarios para almacenar los datos de cada sección
    difusor = {}
    compresor = {}
    cam_comb = {}
    turbina = {}
    tobera = {}

    # Inicializa una variable para rastrear la sección actual
    seccion_actual = None

    # Recorre cada línea del archivo
    for linea in lineas:
        # Elimina espacios en blanco al inicio y al final de la línea
        linea = linea.strip()

        # Comprueba si la línea indica una nueva sección de datos
        if linea.startswith('#DATOS'):
            # Obtiene el nombre de la sección eliminando '#DATOS '
            seccion = linea.replace('#DATOS ', '')
            # Actualiza la sección actual
            seccion_actual = seccion.lower().replace(' ', '_')
            # Continúa con la siguiente línea
            continue
        
        # Si estamos en una sección válida y no es una línea de comentario
        if seccion_actual and not linea.startswith('#'):
            # Divide la línea en nombre de variable y valor utilizando el signo "=" como separador
            nombre, valor = linea.split('=')
            # Almacena el nombre de la variable y su valor en el diccionario correspondiente a la sección actual
            if seccion_actual == 'difusor':
                difusor[nombre.strip()] = valor.strip()
            elif seccion_actual == 'compresor':
                compresor[nombre.strip()] = valor.strip()
            elif seccion_actual == 'camara_comb':
                cam_comb[nombre.strip()] = valor.strip()
            elif seccion_actual == 'turbina':
                turbina[nombre.strip()] = valor.strip()
            elif seccion_actual == 'tobera':
                tobera[nombre.strip()] = valor.strip()
        
        #Se crea una lista de dicccionarios para acceder a ellos posteriormente
        variables = [difusor, compresor, cam_comb, turbina, tobera]
    # Devuelve un diccionario que contiene los datos de todas las secciones
    return variables

def abrir_txt_turb(): #Abre el archivo de datos de la turbina
    # Ruta a la carpeta donde está el archivo
    carpeta = os.path.join('..', 'datos')
    nombre_archivo = 'turbine_data.txt'

    # Ruta completa al archivo
    archivo = os.path.join(carpeta, nombre_archivo)
    return archivo

def abrir_txt_resultados(): #Abre el archivo de resultados
    # Ruta a la carpeta donde está el archivo
    carpeta = os.path.join('..', 'resultados', 'rpm', 'tobera_conv_div')
    nombre_archivo = 'resultados_supersónico.txt'

    # Ruta completa al archivo
    archivo = os.path.join(carpeta, nombre_archivo)
    return archivo

def guardar_resultados_txt(archivo, num_prueba, variable, valor, resultados):
    with open(archivo, 'a') as f:
        f.write(f"Prueba número {num_prueba} variando la {variable}\n")
        f.write(f"{variable} = {valor}\n")
        f.write(f"Velocidad de vuelo = {resultados[0]}\n")
        f.write(f"Potencia compresor = {resultados[1]}\n")
        f.write(f"Potencia turbina = {resultados[2]}\n")
        f.write(f"Empuje = {resultados[3]}\n")
        f.write(f"Masa de combustible = {resultados[4]}\n")
        f.write(f"Impulso específico = {resultados[5]}\n")
        f.write(f"TSFC = {resultados[6]}\n")
        f.write(f"Tobera = {resultados[7]}\n")
        if variable == 'altura':
            f.write(f"rpm = {resultados[8]}\n")
        elif variable == 'RPM':
            f.write(f"altura = {resultados[9]}\n")
        else:
            raise('Introduzca una variable válida: altura o RPM.')
        f.write("\n")
main

def secant_method(func, x0, x1, max_iter=1000):
    """
    Encuentra una raíz de la función 'func' utilizando el método de la secante.
    
    Args:
        func: La función para la cual se busca la raíz.
        x0: Primer punto inicial.
        x1: Segundo punto inicial.
        tol: Tolerancia para el criterio de parada (por defecto es 1e-6).
        max_iter: Número máximo de iteraciones permitidas (por defecto es 1000).
        
    Returns:
        La aproximación de la raíz encontrada.
    """
    iter_count = 0
    while iter_count < max_iter:
        x2 = x1 - func(x1)[0] * ((x1 - x0) / (func(x1)[0] - func(x0)[0]))
        print(x2, 'x2 esta es la temeperatura siguiente')
        dif = func(x2)[0]
        print(dif,'diferencia de potencias')
        if dif < 0.00001:
            return func(x2)
        x0, x1 = x1, x2
        iter_count += 1
    print("El método de la secante no convergió después de", max_iter, "iteraciones.")
    return None

def main():
    # Llamada a la función 'difusor' pasando como argumento el resultado de 'abrir_txt()', que se espera que devuelva algún tipo de datos.
    resultados_difusor = difusor(abrir_txt())
    
    # Llamada a la función 'solver_compresor' con los resultados del difusor y los datos de 'abrir_txt()'.
    resultados_compresor = solver_compresor(abrir_txt(), resultados_difusor)
    
    # Obtención de la potencia del compresor y la eficiencia mecánica.
    W_comp = float(resultados_compresor['Potencia'])
    eta_m = float(abrir_txt()[3]['eta_m'])
    
    # Definición de una función interna 'Potencias' que toma un parámetro 'T_out_cc'.
    def Potencias(T_out_cc):
        # Llamada a la función 'cam_comb' con resultados del compresor y 'T_out_cc', obteniendo los resultados de la cámara de combustión y el combustible.
        resultados_camcomb, comp_cc = cam_comb(abrir_txt(), resultados_compresor, T_out_cc)
        
        # Llamada a la función 'main_turbina' con varios parámetros, obteniendo los resultados de la turbina.
        resultados_turbina = main_turbina(abrir_txt()[1], abrir_txt()[2], resultados_camcomb, abrir_txt_turb())
        
        # Obtención de la potencia de la turbina.
        W_turb = float(resultados_turbina['Potencia'])
        
        # Cálculo de la diferencia entre la potencia del compresor y la potencia de la turbina ajustada por la eficiencia mecánica.
        dif = W_comp - W_turb * eta_m
        
        # Retorna múltiples valores calculados.
        return dif, resultados_turbina, comp_cc, resultados_camcomb, W_turb
    
    # Uso del método de la secante para encontrar el valor de 'T_out_cc' que hace que la función 'Potencias' sea cero, dentro del rango 900 a 1200.
    dif, resultados_turbina, comp_cc, resultados_camcomb, W_turb = secant_method(Potencias, 900, 1200 )
    
    # Llamada a la función 'main_tobera' con los resultados de las etapas anteriores.
    resultados_tobera = main_tobera(abrir_txt(), resultados_difusor, resultados_turbina, comp_cc)
    
    # Extracción de varios resultados de la tobera.
    tobera = str(resultados_tobera['tobera'])
    m_dot = float(resultados_tobera['m_dot'])
    P_out = float(resultados_tobera['P_out'])
    c_out = float(resultados_tobera['c_out'])
    P_e = float(resultados_tobera['P_e'])
    A_out = float(resultados_tobera['A_out'])
    c_in = float(abrir_txt()[0]['v_vuelo'])
    H = float(abrir_txt()[0]['H_vuelo'])
    rpm = float(abrir_txt()[1]['w_rpm'])
    
    # Cálculo del empuje.
    T = m_dot * (c_out - c_in) + A_out * (P_out - P_e)
    
    # Obtención del flujo de masa de combustible.
    m_cble = float(resultados_camcomb['m_cble'])
    
    # Cálculo del consumo específico de combustible.
    TSFC = m_cble / T
    
    # Cálculo del impulso específico.
    I_s = 1 / (9.81 * TSFC)
    
    return c_in, W_comp, W_turb, T, m_cble, I_s, TSFC, tobera, rpm, H



if __name__ == '__main__':
    resultados = main()
    print(resultados)
    guardar_resultados_txt(abrir_txt_resultados(), 0, 'RPM', abrir_txt()[1]['w_rpm'], resultados)


