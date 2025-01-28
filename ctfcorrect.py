#! /usr/bin/env python

'''
test for ctf correction
'''

import mrcfile
import os
import cupy as cp
from matplotlib import pyplot


def create_astigmatic_ctf(size, defocus_u, defocus_v, astig_angle, voltage, cs, amplitude_contrast, pixel_size):
    # Create frequency grid
    freq_x = fftpack.fftfreq(size, d=pixel_size)
    freq_y = fftpack.fftfreq(size, d=pixel_size)
    fx, fy = np.meshgrid(freq_x, freq_y)
    
    # Convert astigmatism angle to radians
    angle_rad = np.deg2rad(astig_angle)
    
    # Calculate spatial frequency
    k2 = fx**2 + fy**2
    k = np.sqrt(k2)
    
    # Calculate the azimuthal angle for each pixel
    # Note: We use arctan2 to get the correct quadrant
    az = np.arctan2(fy, fx)
    
    # Calculate effective defocus for each pixel considering astigmatism
    defocus_average = (defocus_u + defocus_v) / 2
    defocus_diff = (defocus_u - defocus_v) / 2
    defocus_z = defocus_average + defocus_diff * np.cos(2 * (az - angle_rad))
    
    # Calculate wavelength (in nm)
    wavelength = 12.26 / np.sqrt(voltage * (1 + voltage * 0.978466e-6))
    
    # Calculate phase shift including astigmatism
    chi = np.pi * wavelength * k2 * (defocus_z - 0.5 * cs * wavelength**2 * k2)
    
    # Create CTF with astigmatism
    ctf = -np.sqrt(1 - amplitude_contrast**2) * np.sin(chi) - amplitude_contrast * np.cos(chi)
    
    return ctf

