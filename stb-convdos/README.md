# stb_convdos.x - DOS Convolution


## Description
This program applies a Gaussian convolution to Density of States (DOS) data provided in an input file. The goal is to smooth the data and reduce numerical fluctuations, making it easier to interpret electronic structure results.

Gaussian convolution is commonly used in computational physics and materials science to remove noise from DOS calculations, particularly those obtained from Density Functional Theory (DFT) simulations.

## Installation
The script is written in Python 3 and depends on the following libraries:

- `numpy` (for numerical computations)
- `matplotlib` (for plotting results)

To install the required dependencies, run:
```bash
pip install numpy matplotlib
```

## Usage
Run the script with the following parameters:
```bash
stb_convdos.x --input-file dos.dat --size 100 --sigma 1.0 --outfile dos_conv.dat
```

### Parameters:
- `--input-file` (required): Name of the input file containing DOS data.
- `--size` (required): Size of the Gaussian mask (number of points used in smoothing).
- `--sigma` (required): Standard deviation of the Gaussian function, controlling the degree of smoothing.
- `--outfile` (required): Name of the output file for the filtered data.

### Example usage:
```bash
stb_convdos.x --input-file dos.dat --size 50 --sigma 0.5 --outfile dos_smooth.dat
```

## Input Format
The input file should be a text file with two or more columns:
- The first column represents energy values.
- The second (and additional) columns contain the corresponding DOS values.
- The script processes each DOS column independently.

Example of a valid input file (`dos.dat`):
```txt
-5.0  0.002
-4.5  0.005
-4.0  0.009
...
5.0   0.001
```

## Output
- A text file containing the smoothed DOS data, maintaining the same structure as the input file.
- A graph saved in `.png` format displaying both the original and smoothed DOS curves for easy comparison.

## Example Output File (`dos_smooth.dat`):
```txt
-5.0  0.001
-4.5  0.004
-4.0  0.008
...
5.0   0.0009
```

## Applications
This tool is particularly useful for:
- Post-processing DOS data obtained from **DFT calculations**.
- Improving the visualization of electronic structure features.
- Removing numerical artifacts and oscillations in the DOS curves.

## Author
**Dr. Carlos M. O. Bastos**  
University of Brasilia (UnB), Brazil

## License
This software is distributed with stb-suite

## References
For more details on DOS calculations and Gaussian convolution, consider:
- Richard M. Martin, *Electronic Structure: Basic Theory and Practical Methods* (Cambridge University Press, 2004).
- D. S. Sholl and J. A. Steckel, *Density Functional Theory: A Practical Introduction* (Wiley, 2009).

## Contributing
Contributions are welcome! If you'd like to improve this project, please submit a pull request or open an issue.

## Contact
For questions or suggestions, please contact **Dr. Carlos M. O. Bastos** at [carlos.bastos@unb.br](mailto:carlos.bastos@unb.br).

