#!/usr/bin/python
# Script for plotting spectra based on a .dat file containing a list of vaules and heights.
# Tested with Python 3.8.16

import numpy as np

from argparse import ArgumentParser
# from math import pi

def read_command_line_arguments():
    # Uses argparse library to read and assign 
    parser = ArgumentParser(description="This script processes a set of frequencies and intensities to output an IR spectrum with broadened peaks.")
    parser.add_argument('input_filename', type=str, help="File containing a set of IR frequencies (in cm^-1) and corresponding intensities")
    parser.add_argument('output_filename', type=str, help="File where the spectrum is written")
    parser.add_argument("-lw", "--linewidth", type=float, default=5.,help="Gaussian/Lorentzian broadening parameter. Defalut: 5.")
    parser.add_argument("-s", "--shift", type=float, default=1.,help="Scaling factor applied to the spectrum. Default: 1.")
    parser.add_argument("-r", "--range", type=list, default=[0,3500],help="Range of values over which the output spectrum is plotted. Default: [0,3500]")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-G", "--gaussian", action="store_true", help="Use Gaussian peak broadening.")
    group.add_argument("-L", "--lorentzian", action="store_true", help="Use Lorentzian peak broadening. This is the default setting.")
    arguments = parser.parse_args()
    return arguments

def read_file(filename) -> [list, list]:
    #Parses the file containing two columns of numerical data
    
    #Parameters:
    #    filename (str) : Data file (with or without header) containing two columns

    #Returns:
    #    x (list) : Values that denote the position of the peak on x-axis
    #    y (list) : Values that denote the height of the corresponding peak

    with open(filename, 'r') as file:
        x = []
        y = []
        for line in file:
            p = line.split()
            x.append(float(p[0]))
            y.append(float(p[1]))

    return x, y

def write_file(filename,array):
   #Writes out the resulting spectrum to a .dat file

   #Parameters:
   #    filename (str) : Output datafile containing spectrum
   #    array (np.Array) : The result of the spectrum generation, with intensity values at even intervals.
   
   np.savetxt(filename,array,fmt='%.2f %0.6f',delimiter='\t')

def generate_spectrum(freq,intensity,pw):
    # Generates an IR spectrum for a set of frequencies and associated intensities 
    
    # Parameters:
    #    freq (list) : set of IR frequencies corresponding to centre of peak
    #    inten (float) : corresponding set of IR intensities for the IR peaks
    #    pw (float) : peak width parameter (sigma if gaussian, tau if lorenzian)

    # Returns:
    #    spectrum (np.array) : a set of absorption intensities associated with x-axis values in spectrum_range
    
    if shift != 1.:
        freq = list(map(lambda x: x * shift, freq))

    spectrum_range = np.arange(plot_range[0],plot_range[1],1.) # Resolution of output spectrum is set to 1. here, and can be altered.
    spectrum_range = np.append(spectrum_range,freq)
    spectrum_range.sort()
    spectrum = np.zeros(len(spectrum_range))

    if lorentz and not gauss:
        broadening_kernel = lorentzian
    elif gauss and not lorentz:
        broadening_kernel = gaussian
    else:
        raise ValueError("Please specify either Gaussian or Lorentzian broadening.")
    
    for freq,intensity in zip(freq,intensity):
        min_peak = freq - 20.*pw # The 20 can be changed if the spectrum starts 
        max_peak = freq + 20.*pw
        one_peak = [0.0 if i < min_peak else 0.0 if i > max_peak else intensity*broadening_kernel(i,freq,pw) for i in spectrum_range]
        spectrum += np.asarray(one_peak)

    return np.column_stack((spectrum_range,spectrum))

def lorentzian(x, x0, tau : float) -> float :
    # Generates value of a Lorentzian function centered at x0 with linewidth tau at coordinate x.
    
    # Parameters:
    #    x (list) : x value at which the spectrum is generated
    #    x0 (float) : centre of the lorenzian function
    #    tau (float) : parameter of lorenzian function broadening

    return (tau**2)/((x-x0)**2+(tau**2))

def gaussian(x, x0, sigma : float) -> float :
    # Generates value of a Gaussian function centered at x0 with sigma at coordinate x.
    
    # Parameters:
    #    x (list) : x value at which the spectrum is generated
    #    x0 (float) : centre of the gaussian function
    #    sigma (float) : parameter of gaussian function broadening

    return np.exp(-np.power((x - x0)/sigma, 2.)/2.)

if __name__ == "__main__":
    args = read_command_line_arguments()

    infilename = args.input_filename
    outfilename = args.output_filename
    linewidth = args.linewidth
    shift = args.shift
    gauss = args.gaussian
    lorentz = args.lorentzian
    plot_range = args.range

    xval, yval = read_file(infilename)
    spectrum = generate_spectrum(xval,yval,linewidth)
    write_file(outfilename,spectrum)



