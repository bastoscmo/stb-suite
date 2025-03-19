#!/usr/bin/env python3
import os
import subprocess

# Define input file and format
input_file = "test/test.in"
input_format = "poscar"  # Assuming input is in POSCAR format
output_formats = ["xyz", "poscar", "fdf", "dftb", "xsf", "fhi"]

# Loop through all output formats
for out_format in output_formats:
    output_file = f"test/convert/test.{out_format}"

    # Build the command
    command = [
        "python3", "stb-translate.x",
        "-if", input_format,
        "-i", input_file,
        "-of", out_format,
        "-o", output_file
    ]

    # Execute the command
    print(f"Converting {input_file} to {output_file}...")
    result = subprocess.run(command, capture_output=True, text=True)

    # Print results
    if result.returncode == 0:
        print(f"✅ Successfully converted to {out_format}")
    else:
        print(f"❌ Error converting to {out_format}: {result.stderr}")

print("✅ All conversions completed!")
