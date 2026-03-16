import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.fftpack import fft2, ifft2

# --- 1. Parámetros Técnicos (Ajustados para visibilidad) ---
N = 128                 # Resolución equilibrada
L = 2 * np.pi
dt = 20        # Paso un poco más largo para ver movimiento
visc = 1e-6             # Viscosidad ultra-baja para evitar que se "borre"
p_nu = 1                # Operador de difusión estándar (nabla^2)

dx = L / N
x = np.linspace(0, L, N, endpoint=False)
y = np.linspace(0, L, N, endpoint=False)
X, Y = np.meshgrid(x, y)

# Frecuencias
k = np.fft.fftfreq(N, dx) * 2 * np.pi
l = np.fft.fftfreq(N, dx) * 2 * np.pi
K, L_k = np.meshgrid(k, l)
K_sq = K**2 + L_k**2
inv_K_sq = -1.0 / K_sq
inv_K_sq[0, 0] = 0.0  # Evitar división por cero

# --- 2. Inicialización de Energía ---
# En lugar de ruido puro, usamos un espectro de potencia que favorece remolinos
q_hat = np.complex128(np.random.randn(N, N) + 1j * np.random.randn(N, N))
q_hat *= np.exp(-K_sq / (2 * (10.0**2))) # Filtro para crear manchas visibles
q_hat[0, 0] = 0.0

fig, ax = plt.subplots(figsize=(8, 8), facecolor='black')
# Usamos 'RdBu_r' para que el rojo sea positivo y azul negativo
im = ax.imshow(np.zeros((N, N)), cmap='RdBu_r', origin='lower', extent=[0, L, 0, L])
ax.axis('off')
txt = ax.text(0.05, 0.95, "", color='white', transform=ax.transAxes, fontfamily='monospace')

def get_dqdt(q_hat_loc):
    # Inversión de PV
    psi_hat = q_hat_loc * inv_K_sq
    
    # Velocidades y Gradientes
    u = -ifft2(1j * L_k * psi_hat).real
    v =  ifft2(1j * K * psi_hat).real
    dqdx = ifft2(1j * K * q_hat_loc).real
    dqdy = ifft2(1j * L_k * q_hat_loc).real
    
    # Advección y Disipación
    Jacobian = u * dqdx + v * dqdy
    return -fft2(Jacobian) - visc * K_sq * q_hat_loc

def update(frame):
    global q_hat
    
    # Integración RK4
    for _ in range(15):
        k1 = get_dqdt(q_hat)
        k2 = get_dqdt(q_hat + 0.5 * dt * k1)
        k3 = get_dqdt(q_hat + 0.5 * dt * k2)
        k4 = get_dqdt(q_hat + dt * k3)
        q_hat += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    q_real = ifft2(q_hat).real
    im.set_data(q_real)
    
    # Normalización dinámica: esto evita la pantalla negra
    vmax = np.max(np.abs(q_real))
    im.set_clim(-vmax, vmax)
    
    txt.set_text(f"Paso: {frame} | Energía Máx: {vmax:.4f}")
    return [im, txt]

ani = animation.FuncAnimation(fig, update, frames=400, interval=20, blit=True)
plt.show()
  
