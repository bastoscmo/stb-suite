#!/usr/bin/env python

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################

VERSION = "1.8.0"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import numpy as np
import argparse
import matplotlib.pyplot as plt


# ANSI colors for terminal
COLORS = {
    'reset': '\033[0m',
    'cyan': '\033[96m',
    'blue': '\033[94m',
    'green': '\033[92m',
    'yellow': '\033[93m',
    'red': '\033[91m',
    'bold': '\033[1m',
    'underline': '\033[4m'
}

def color_text(text: str, color: str) -> str:
    """Returns text formatted with ANSI color"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Displays the stylized STB-SUITE intro"""
    os.system('cls' if os.name == 'nt' else 'clear')
    
    logo = color_text(r"""
.----------------.  .----------------.  .----------------.
| .--------------. || .--------------. || .--------------. |
| |    _______   | || |  _________   | || |   ______     | |
| |   /  ___  |  | || | |  _   _  |  | || |  |_   _ \    | |
| |  |  (__ \_|  | || | |_/ | | \_|  | || |    | |_) |   | |
| |   '.___`-.   | || |     | |      | || |    |  __'.   | |
| |  |`\____) |  | || |    _| |_     | || |   _| |__) |  | |
| |  |_______.'  | || |   |_____|    | || |  |_______/   | |
| |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'
 """, 'cyan')

    description = [
        "Siesta ToolBox Suite",
        "A comprehensive toolkit for SIESTA DFT simulations",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.2)
    print("="*60 + "\n")
    return



def gaussian_kernel_1d(size, sigma):
    """Creates a 1D Gaussian kernel."""
    kernel_1d = np.linspace(-(size // 2), size // 2, size)
    kernel_1d = np.exp(-(kernel_1d**2) / (2 * sigma**2))
    kernel_1d /= np.sum(kernel_1d)  # Normalize the kernel
    return kernel_1d

def convolve1d(data, kernel):
    """Performs 1D convolution of a vector with a kernel."""
    kernel_size = kernel.shape[0]
    data_size = data.shape[0]
    pad_size = kernel_size // 2
    padded_data = np.pad(data, (pad_size, pad_size), mode='constant')
    output = np.zeros_like(data)
    for i in range(data_size):
        region = padded_data[i:i + kernel_size]
        output[i] = np.sum(region * kernel)
    return output

# --- CORRECTED PLOT FUNCTION ---
def plot(energy_data, original_dos_data, filtered_dos_data, col_index):
    """
    Documentation:
    Plots the original and filtered data for a specific DOS column.
    
    Arguments:
    energy_data (np.array): Vector with energy values (X-axis).
    original_dos_data (np.array): Vector with original DOS data (Y-axis of 1st plot).
    filtered_dos_data (np.array): Vector with filtered DOS data (Y-axis of 2nd plot).
    col_index (int): The index of the column being plotted (for the title).
    """
    plt.figure(figsize=(10, 5))
    
    # Subplot 1: Original Data
    plt.subplot(1, 2, 1)
    # Documentation: The title now shows which column (by index) is being displayed.
    plt.title(f"Original")
    # Documentation: Plots energy vs. the specific original DOS column.
    plt.plot(energy_data, original_dos_data)
    plt.xlabel("Energy")
    plt.ylabel("DOS")

    # Subplot 2: Filtered Data
    plt.subplot(1, 2, 2)
    # Documentation: Corresponding title for the filtered data.
    plt.title(f"Filtered")
    # Documentation: Plots energy vs. the specific filtered DOS column.
    plt.plot(energy_data, filtered_dos_data)
    plt.xlabel("Energy")
    plt.ylabel("DOS")
    
    # Documentation: Shows the plot.
    plt.show()

# --- CORRECTED FILTER_DATA FUNCTION ---
def filter_data(inp_file, kernel):
    """
    Documentation:
    Applies convolution to all data columns of the input file
    and plots each step.
    
    Arguments:
    inp_file (np.array): The complete array read from the file.
    kernel (np.array): The Gaussian kernel for convolution.
    
    Returns:
    np.array: A new array containing the energy column and all
              filtered DOS columns.
    """
    
    # Documentation: Separates the energy column (assuming it's column 0).
    energy_data = inp_file[:, 0]
    
    # Documentation: Initializes the filtered data array ONLY with the energy column.
    # We will build this array column by column.
    # We use reshape to ensure it is a column vector
    filtered_data_final = energy_data.reshape(-1, 1)
    
    # Documentation: Iterates over the data columns, starting from column 1
    # (since column 0 is energy).
    # We use inp_file.shape[1] to get the total number of columns.
    for i in range(1, inp_file.shape[1]):
        
        # Documentation: Selects the current original DOS column (e.g., column 1, then 2, etc.)
        original_dos_data = inp_file[:, i]
        
        # Documentation: Applies convolution ONLY to this specific column.
        filtered_dos_data = convolve1d(original_dos_data, kernel)
        
        # Documentation: Adds (stacks) the new filtered column to our final array.
        # We use reshape to ensure the vector is a column before stacking.
        filtered_data_final = np.column_stack((filtered_data_final, filtered_dos_data.reshape(-1, 1)))
        
        # Documentation: Calls the UPDATED plot function.
        # We pass (energy, current_original_dos, current_filtered_dos, column_index)
        plot(energy_data, original_dos_data, filtered_dos_data, i)
        
    # Documentation: Returns the complete array with the energy and all filtered columns.
    return filtered_data_final
# --- END OF CORRECTIONS ---

def main():
    parser = argparse.ArgumentParser(description="Apply the Gaussian Convolution in DOS.")
    parser.add_argument("--file", dest="input_file",required=True, help="Input file with  DOS.")
    parser.add_argument("--size", type=int, required=True, help="Size of Gaussian mask.")
    parser.add_argument("--sigma", type=float, required=True, help="Standard deviation of the Gaussian function.")
    parser.add_argument("--out", required=True,dest="outfile", help="Output file with filtered data.")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")
    args = parser.parse_args()

    if args.intro == True:
        show_intro()
    print("\n" + color_text("DOS Convolution Tool:", 'bold'))
    print("-"*60)

    kernel = gaussian_kernel_1d(args.size, args.sigma)
    inp_file = np.loadtxt(args.input_file)
    print("\n[INFO] Read File") 
    print("[INFO] Applied DOS convolution and plotting...") 
    print("[WARNING] \n") 
    filtered_data = filter_data(inp_file, kernel)
    np.savetxt(args.outfile, filtered_data, fmt='%.6f', header="Energy DOS_filtered") 
    print(f"\n[OK] Filtered data write in {args.outfile}")
    print("[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("We’ve convoluted everything… including the soul of the electron.\n\n", 'bold'))

if __name__ == "__main__":
    main()
