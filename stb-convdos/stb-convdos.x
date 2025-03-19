#!/usr/bin/env python3

#################################################
# Siesta Tool Box - Convolution of the DOS      #
# version 1.0.0                                 #
# University of Brasilia (UnB) - Brazil         #
# version 1.0.0 2025/03/06                      #
# Dr. Carlos M. O. Bastos                       #
#################################################

__version__ = "1.0.0"

import numpy as np
import argparse
import matplotlib.pyplot as plt

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
    parser.add_argument("--input-file", required=True, help="Input file with  DOS.")
    parser.add_argument("--size", type=int, required=True, help="Size of Gaussian mask.")
    parser.add_argument("--sigma", type=float, required=True, help="Standard deviation of the Gaussian function.")
    parser.add_argument("--outfile", required=True, help="Output file with filtered data.")
    args = parser.parse_args()

    kernel = gaussian_kernel_1d(args.size, args.sigma)
    inp_file = np.loadtxt(args.input_file)
    filtered_data = filter_data(inp_file, kernel)

    np.savetxt(args.outfile, filtered_data, fmt='%.6f', header="Energy DOS_filtered")
    print(f"Processo concluído. Dados filtrados salvos em {args.outfile}")

if __name__ == "__main__":
    main()
