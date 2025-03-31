import numpy as np
import matplotlib.pyplot as plt

# Función para leer el archivo de datos
def leer_datos(filename):
    # Leer los datos desde el archivo
    data = np.loadtxt(filename)
    
    # Asumiendo que las columnas son: 
    # [Temperatura, Magnetización, Magnetización^2, Error, mc*tau, Correlación]
    T = data[:, 0]  # Temperatura
    rm = data[:, 1]  # Magnetización promedio
    rm2 = data[:, 2]  # Magnetización^2 promedio
    error = data[:, 3]  # Error
    mc_tau = data[:, 4]  # mc*tau
    c = data[:, 5]  # Correlación
    
    return T, rm, rm2, error, mc_tau, c

# Función para graficar los resultados
def graficar_resultados(T, rm, error, c):
    # Crear una figura para los gráficos
    fig, ax1 = plt.subplots(figsize=(8, 6))

    # Graficar la magnetización promedio
    ax1.set_xlabel('Temperatura (T)')
    ax1.set_ylabel('Magnetización Promedio', color='tab:blue')
    ax1.plot(T, rm, color='tab:blue', label='Magnetización', linestyle='-', marker='o')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    
    # Crear un segundo eje y para la correlación y error
    ax2 = ax1.twinx()  
    ax2.set_ylabel('Correlación', color='tab:red')
    ax2.plot(T, c, color='tab:red', label='Correlación', linestyle='--')
    ax2.tick_params(axis='y', labelcolor='tab:red')
    
    # Graficar el error
    ax2.plot(T, error, color='tab:green', label='Error', linestyle='-.')
    
    # Título y leyendas
    plt.title('Magnetización y Correlación en función de la Temperatura')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    
    # Mostrar el gráfico
    plt.show()

# Función principal
def main():
    """
    A ver qué obtengo con el primer fichero
    """
    # Ruta al archivo de salida generado por Fortran
    filename = 'fort.66'
    
    # Leer los datos
    T, rm, rm2, error, mc_tau, c = leer_datos(filename)
    
    # Graficar los resultados
    graficar_resultados(T, rm, error, c)

    """
    A ver qué obtengo con el segundo
    """
    # Leer los datos de fort.88 (magnetización)
    filename = 'fort.88'

    # Cargar los datos de magnetización (se asume que son valores de punto flotante por línea)
    magnetization_data = np.loadtxt(filename)

    # Graficar los datos
    plt.figure(figsize=(10, 6))
    plt.plot(magnetization_data, marker='o', linestyle='-', color='r', label='Magnetización')
    plt.title("Evolución de la magnetización durante la simulación")
    plt.xlabel("Paso de simulación")
    plt.ylabel("Magnetización")
    plt.legend()
    plt.grid(True)
    plt.show()

    archivo = "fort.88"
    magnetizacion = np.loadtxt(archivo)

    # Graficar la magnetización en función de la iteración
    plt.figure(figsize=(8, 5))
    plt.plot(magnetizacion, marker='.', linestyle='-', alpha=0.7)
    plt.xlabel("Iteraciones de medición")
    plt.ylabel("Magnetización instantánea |M|/N")
    plt.title("Evolución de la magnetización en el tiempo")
    plt.grid()
    plt.show()

# Llamada a la función principal
if __name__ == '__main__':
    main()
