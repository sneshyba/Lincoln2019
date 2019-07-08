
# Does 2d Fourier transform (to/from) 

import numpy as np
def FT(grid,x,y):


    # Fourier transform it
    grid_FT = np.fft.fft2(grid)
    grid_FTshift = np.fft.fftshift(grid_FT)

    # And get the frequencies
    Ny, Nx = np.shape(grid)
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    kx = np.fft.fftfreq(Nx)*2*np.pi/dx
    ky = np.fft.fftfreq(Ny)*2*np.pi/dy
    kxshift = np.fft.fftshift(kx)
    kyshift = np.fft.fftshift(ky)
    return grid_FTshift, kxshift, kyshift


def IFT(sollast_FTshift_filtered):

    # Un-shift it
    sollast_FT_filtered = np.fft.ifftshift(sollast_FTshift_filtered)

    # Inverse Fourier transform
    sollast_FT_filtered_IFT = np.fft.ifft2(sollast_FT_filtered)

    return sollast_FT_filtered_IFT
