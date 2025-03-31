import numpy as np
import matplotlib.pyplot as plt

# Función para leer el archivo de datos
def leer_datos(filename):
    data = np.loadtxt(filename)
    T, rm, rm2, error, mc_tau, c = data.T  # Transponer para asignar correctamente
    return T, rm, rm2, error, mc_tau, c

# Cargar los datos del primer archivo (fort.66)
filename1 = 'fort.66'
temperatura, magnetizacion, magnetizacion2, error, mc_tau, correlacion = leer_datos(filename1)

# Cargar los datos del segundo archivo (fort.88)
filename2 = 'fort.88'
magnetization_data = np.loadtxt(filename2)  # Asume que son valores de punto flotante por línea

# Parámetros adicionales
L_values = [20, 40, 80, 120, 160, 200]  # Tamaños de red usados en las figuras
Tc = 2.269  # Temperatura crítica del modelo de Ising 2D
J_k = 1.0  # Factor J/k para la ecuación de Onsager

# Solución exacta de Onsager para la magnetización
def onsager_magnetization(T):
    if T < Tc:
        return (1 - (np.sinh(2 * J_k / T) ** -4)) ** (1/8)
    else:
        return 0

osanger_T = np.linspace(0, 4, 500)
onsager_m = np.array([onsager_magnetization(T) for T in osanger_T])

# Figura 8: Magnetización vs Temperatura
plt.figure(figsize=(8,6))
plt.errorbar(temperatura, magnetizacion, yerr=error, fmt='o', label="Simulación")
plt.plot(osanger_T, onsager_m, '-', color='black', label="Onsager Exacto")
plt.axvline(Tc, linestyle='--', color='red', label="$T_c$ exacto")
plt.xlabel("Temperatura T")
plt.ylabel("Magnetización ⟨M⟩")
plt.legend()
plt.title("Figura 8: Magnetización vs Temperatura")
plt.grid()
plt.show()

# Figura 9: Susceptibilidad vs Temperatura
susceptibilidad = (magnetizacion2 - magnetizacion**2) / temperatura
plt.figure(figsize=(8,6))
plt.plot(temperatura, susceptibilidad, 'o-', label="$\chi_T$ (Simulación)")
plt.axvline(Tc, linestyle='--', color='red', label="$T_c$ exacto")
plt.xlabel("Temperatura T")
plt.ylabel("Susceptibilidad $\chi_T$")
plt.legend()
plt.title("Figura 9: Susceptibilidad vs Temperatura")
plt.grid()
plt.show()

# Figura 10: Magnetización reescalada
beta_nu = 1/8 / (1.0)  # Ajusta según corresponda
plt.figure(figsize=(8,6))
plt.plot((1-temperatura/Tc), magnetizacion * L_values[2]**beta_nu, 'o-', label="Reescalado")
plt.xlabel("$(1 - T/T_c)$")
plt.ylabel("$M_0 L^{\beta/\nu}$")
plt.legend()
plt.title("Figura 10: Magnetización reescalada")
plt.grid()
plt.show()

# Figura 11: Susceptibilidad reescalada
plt.figure(figsize=(8,6))
plt.plot((1-temperatura/Tc), susceptibilidad * L_values[2]**(7/4), 'o-', label="Reescalado")
plt.xlabel("$(1 - T/T_c)$")
plt.ylabel("$\chi L^{\gamma/\nu}$")
plt.legend()
plt.title("Figura 11: Susceptibilidad reescalada")
plt.grid()
plt.show()

# Figura 12: Correlación vs Temperatura (si aplica)
plt.figure(figsize=(8,6))
plt.plot(temperatura, correlacion, 'o-', label="Correlación")
plt.xlabel("Temperatura T")
plt.ylabel("Correlación")
plt.legend()
plt.title("Figura 12: Correlación vs Temperatura")
plt.grid()
plt.show()