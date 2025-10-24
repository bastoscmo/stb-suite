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


# Cores ANSI para terminal
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
    """Retorna texto formatado com cor ANSI"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"

def show_intro() -> None:
    """Exibe a introdução estilizada da STB-SUITE"""
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
    """Cria uma máscara gaussiana 1D."""
    kernel_1d = np.linspace(-(size // 2), size // 2, size)
    kernel_1d = np.exp(-(kernel_1d**2) / (2 * sigma**2))
    kernel_1d /= np.sum(kernel_1d)  # Normaliza a máscara
    return kernel_1d

def convolve1d(data, kernel):
    """Realiza a convolução 1D de um vetor com um kernel."""
    kernel_size = kernel.shape[0]
    data_size = data.shape[0]
    pad_size = kernel_size // 2
    padded_data = np.pad(data, (pad_size, pad_size), mode='constant')
    output = np.zeros_like(data)
    for i in range(data_size):
        region = padded_data[i:i + kernel_size]
        output[i] = np.sum(region * kernel)
    return output

def plot(inp_file, filtered_data):
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.title("Original")
    plt.plot(inp_file[:, 0], inp_file[:, 1])
    plt.xlabel("Energy")
    plt.ylabel("DOS")

    plt.subplot(1, 2, 2)
    plt.title("Filtered")
    plt.plot(filtered_data[:, 0], filtered_data[:,-1])
    plt.xlabel("Energy")
    plt.ylabel("DOS")
    plt.show()

def filter_data(inp_file, kernel):
    filtered_data=inp_file[:,0]
    for i in range(1,len(inp_file.T)):
        filtered_data = np.column_stack((filtered_data, convolve1d(inp_file[:, i], kernel)))
        plot(inp_file,filtered_data)
    return filtered_data

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
