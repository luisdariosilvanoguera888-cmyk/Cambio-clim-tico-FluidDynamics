# 🌊 2D Navier-Stokes Pseudo-Spectral Solver

Este proyecto es un simulador de dinámica de fluidos incompresibles en 2D que utiliza un enfoque **pseudo-espectral** para resolver la ecuación de transporte de vorticidad. Es ideal para estudiar la turbulencia en condiciones ideales y visualizar la formación de estructuras coherentes (vórtices).

![Fluid Simulation Demo](https://github.com/luisdariosilvanoguera888-cmyk/Cambio-clim-tico-FluidDynamics/blob/main/VID-20260316-WA00452-ezgif.com-video-to-gif-converter.gif) ## 🔬 Fundamentos Físicos

El simulador resuelve la evolución de la vorticidad $q$ definida como el rotacional de la velocidad ($\omega = \nabla \times \mathbf{u}$). La dinámica se rige por:

$$\frac{\partial q}{\partial t} + \mathbf{u} \cdot \nabla q = \nu \nabla^2 q$$

### Implementación Técnica:
1.  **Espacio de Fourier:** Las derivadas espaciales se calculan en el dominio de las frecuencias mediante la **FFT (Transformada Rápida de Fourier)**, lo que proporciona una precisión espectral superior a los métodos de diferencias finitas.
2.  **Inversión de Poisson:** Para hallar el campo de velocidades a partir de la vorticidad, resolvemos la ecuación de Poisson para la función de corriente $\psi$:
    $$\hat{\psi} = -\frac{\hat{q}}{k^2 + l^2}$$
3.  **Integración Temporal:** Se utiliza un esquema de **Runge-Kutta de 4to orden (RK4)** para garantizar la estabilidad numérica incluso en pasos de tiempo más largos.

## 🛠️ Características Principales
* **Viscosidad Ultra-Baja:** Optimizado para observar la advección sin una disipación inmediata de la energía.
* **Normalización Dinámica:** Ajuste en tiempo real de la escala de color para visualizar la evolución del sistema a medida que la energía decae.
* **RK4 Integrator:** Implementación manual de integración de alta precisión para sistemas dinámicos.
