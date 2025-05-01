import numpy as np
import matplotlib.pyplot as plt

def leer_configuraciones(nombre_archivo, L):
    """Lee configuraciones asumiendo L valores por línea, L líneas por configuración."""
    configuraciones = []
    valores_por_config = L * L
    valores = []
    linea_num = 0
    
    # Leer todas las líneas y extraer valores
    with open(nombre_archivo, "r") as f:
        for linea in f:
            linea = linea.strip()
            linea_num += 1
            if not linea:
                print(f"Línea {linea_num}: vacía, ignorada")
                continue
            try:
                datos = [int(x) for x in linea.split()]
                print(f"Línea {linea_num}: {len(datos)} valores")
                if len(datos) != L:
                    print(f"Error: esperados {L} valores por línea, encontrados {len(datos)} en línea {linea_num}")
                    return []
                valores.extend(datos)
            except ValueError:
                print(f"Error al leer la línea {linea_num}: {linea}")
                return []
    
    # Verificar el total de valores
    total_valores = len(valores)
    print(f"Total de valores leídos: {total_valores}")
    if total_valores % valores_por_config != 0:
        print(f"Error: el total de valores ({total_valores}) no es múltiplo de {valores_por_config}")
        return []
    
    # Crear configuraciones
    num_configs = total_valores // valores_por_config
    for i in range(num_configs):
        inicio = i * valores_por_config
        fin = inicio + valores_por_config
        config = np.array(valores[inicio:fin]).reshape(L, L)
        configuraciones.append(config)
    
    return configuraciones

def graficar_configuraciones(configuraciones, temperaturas, L):
    n_plots = min(len(configuraciones), len(temperaturas))
    fig, axs = plt.subplots(1, n_plots, figsize=(10, 3))
    if n_plots == 1:
        axs = [axs]
        
    for i, (config, T) in enumerate(zip(configuraciones[:n_plots], temperaturas[:n_plots])):
        axs[i].imshow(config, cmap="gray", interpolation="none")
    #   axs[i].imshow(config, cmap="gray_r", interpolation="none")
        if  T!=2.269: axs[i].set_title(f"$T = {T}$")
        else: axs[i].set_title("$T = T_c$")
        axs[i].axis("off")
    
    plt.tight_layout()
    plt.savefig(fname="lattice.png", dpi=300)
    plt.show()

L = 80
temperaturas = [4.0, 2.269, 2.0, 1.0]

# Leer y graficar
ruta_archivo = "fort.99"  # Ajusta la ruta si es necesario
configuraciones = leer_configuraciones(ruta_archivo, L)
if configuraciones:
    print(f"Se leyeron {len(configuraciones)} configuraciones")
    graficar_configuraciones(configuraciones, temperaturas, L)
else:
    print(f"No se pudieron leer configuraciones válidas de {ruta_archivo}.")