This repository includes code for calculating wavenumber-frequency spectrum and space-time correlation of time series data on a 2-D surface.

The code is implemented in Fortran and relies on [FFTW](https://www.fftw.org/). It supports multi-threaded FFTW for shared-memory machines, enhancing the efficiency of Fourier transformation calculations.

## Definition
### Spectrum
To calculate the wavenumber-frequency spectrum, the time series data $p(x_1,x_3,t)$ (assuming data is zero-centered and homogeneous in $x_1$- and $x_3$-directions) is first Fourier-transformed in space and time,
$$\hat{p}(k_1,k_3,\omega) =\frac{1}{{(2\pi)}^3}\int_{-T/2}^{T/2}\int_{-L_1/2}^{L_1/2}\int_{-L_3/2}^{L_3/2} p(x_1,x_3,t)e^{-(k_1x_1+k_3x_3+\omega t)}dx_1dx_3dt,$$
and then the 2-D wavenumber-frequency spectrum can be calculated as

$$\Phi_{pp}(k_1,k_3,\omega) = \frac{\langle \hat{p}(k_1,k_3,\omega)\hat{p}^*(k_1,k_3,\omega)\rangle}{L_1L_3T},$$

where $\langle\cdot\rangle$ denotes ensemble averaging. From the 2-D spectrum, the 1-D wavenumber-frequency spectrum (for $x_1$-direction) can be calculated as follows:
$$\Phi_{pp}(k_1,\omega) = \int\Phi_{pp}(k_1,k_3,\omega)d k_3.$$

### Correlation

Assuimg data is homogeneous in $x_3$-direction, the space-time correlation $R_{pp}$ in $x_1$-direction is defined as
$$R_{pp}(x_1,\Delta x_1,\Delta t)=\frac{\langle p(x_1,x_3,t)p(x_1+\Delta x_1,x_3,t+\Delta t)\rangle}{\sqrt{\langle p^2(x_1,x_3,t)\rangle}\sqrt{\langle p^2(x_1+\Delta x_1,x_3,t+\Delta t)\rangle}}.$$

Similarly, $R_{pp}$ in $x_3$-direction is defined as
$$R_{pp}(x_1,\Delta x_3,\Delta t)=\frac{\langle p(x_1,x_3,t)p(x_1,x_3+\Delta x_3,t+\Delta t)\rangle}{\sqrt{\langle p^2(x_1,x_3,t)\rangle}\sqrt{\langle p^2(x_1,x_3+\Delta x_3,t+\Delta t)\rangle}}.$$
