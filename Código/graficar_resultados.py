import matplotlib.pyplot as plt
import os

def abrir_txt_resultados(): #Abre el archivo de resultados
    # Ruta a la carpeta donde está el archivo
    carpeta = os.path.join('..', 'resultados','rpm','tobera_convergente')
    nombre_archivo = 'resultados_subsónico.txt'

    # Ruta completa al archivo
    ruta_archivo = os.path.join(carpeta, nombre_archivo)
    return ruta_archivo

def leer_resultados(archivo):
    pruebas = []  # Lista donde se almacenarán las pruebas
    
    # Abrir el archivo en modo lectura
    with open(archivo, 'r') as file:
        lines = file.readlines()  # Leer todas las líneas del archivo y guardarlas en una lista
        
    i = 0
    while i < len(lines):
        if lines[i].startswith("Prueba número"):
            # Extraer el número de prueba, la variable y el valor de la línea actual
            num_prueba = int(lines[i].split()[2])
            variable = lines[i].split()[5]
            valor = float(lines[i+1].split('=')[1].strip())
            
            # Extraer los resultados de las siguientes líneas
            velocidad_vuelo = float(lines[i+2].split('=')[1].strip())
            potencia_compresor = float(lines[i+3].split('=')[1].strip())
            potencia_turbina = float(lines[i+4].split('=')[1].strip())
            empuje = float(lines[i+5].split('=')[1].strip())
            masa_combustible = float(lines[i+6].split('=')[1].strip())
            impulso_especifico = float(lines[i+7].split('=')[1].strip())
            tsfc = float(lines[i+8].split('=')[1].strip())
            tobera = str(lines[i+9].split('=')[1].strip())
            rpm = float(lines[i+10].split('=')[1].strip())

            
            # Crear una lista con los resultados de esta prueba
            resultados = [
                velocidad_vuelo,
                potencia_compresor,
                potencia_turbina,
                empuje,
                masa_combustible,
                impulso_especifico,
                tsfc,
                tobera,
                rpm
            ]
            
            # Agregar la prueba como una tupla a la lista de pruebas
            pruebas.append((num_prueba, variable, valor, resultados))
            
            # Avanzar 12 líneas para saltar al siguiente bloque de datos de prueba
            i += 12
        else:
            # Si la línea actual no comienza con "Prueba número", avanzar a la siguiente línea
            i += 1
    
    return pruebas # Devolver la lista de pruebas procesadas

def graficar_resultados(pruebas, parametro_indice, parametro_nombre, Variable):
    valores = []  # Lista para almacenar los valores de la variable (eje x del gráfico)
    parametros = []  # Lista para almacenar los parámetros específicos (eje y del gráfico)
    
    # Iterar sobre cada prueba en la lista de pruebas
    for prueba in pruebas:
        valor = prueba[2]  # Valor de la variable que se está probando
        parametro = prueba[3][parametro_indice]  # Parámetro específico a graficar
        valores.append(valor)  # Agregar el valor al eje x
        parametros.append(parametro)  # Agregar el parámetro al eje y
    
    # Graficar los datos
    plt.plot(valores, parametros, marker='o')  # Graficar puntos con marcador circular
    plt.xlabel(f'{Variable}')  # Etiqueta del eje x
    plt.ylabel(parametro_nombre)  # Etiqueta del eje y (nombre del parámetro específico)
    #plt.title(f'{parametro_nombre} vs {Variable}')  # Título del gráfico, usando el nombre de la variable probada
    plt.grid(True)  # Mostrar rejilla en el gráfico
    
    # Agregar etiquetas a cada punto
    for i, txt in enumerate(pruebas):
       plt.text(valores[i], parametros[i], f'({int(valores[i])}, {int(parametros[i])})', fontsize=7, ha='right')
   
    plt.show()  # Mostrar el gráfico

    
if __name__ == '__main__':
    archivo = abrir_txt_resultados()
    pruebas = leer_resultados(archivo)
    graficar_resultados(pruebas, 1, 'Potencia Compresor (kW)', 'Velocidad Angular (RPM)')
    graficar_resultados(pruebas, 2, 'Potencia Turbina (kW)', 'Velocidad Angular (RPM)')
    graficar_resultados(pruebas, 3, 'Empuje (N)', 'Velocidad Angular (RPM)')
    #graficar_resultados(pruebas, 4, 'Consumo Combustible (kg/s)', 'Altitud (m)')
    graficar_resultados(pruebas, 5, 'Impulso Específico (s)', 'Velocidad Angular (RPM)')
    #graficar_resultados(pruebas, 6, 'Consumo Específico (kg/N·s)', 'Altitud (m)')
