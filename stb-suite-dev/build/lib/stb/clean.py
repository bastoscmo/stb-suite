#!/usr/bin/env python

#################################################
# Siesta Tool Box - Bands                       #
# version 1.5.1                                 #
# University of Brasilia (UnB) - Brazil         #
# version 1.0.0 2025/02/05                      #
# Dr. Carlos M. O. Bastos                       #
#################################################

VERSION = "1.5.1"

import os
import sys
import warnings
import subprocess
from time import sleep
import argparse
import textwrap
from typing import List, Dict
import numpy as np
import re
import argparse


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

def should_delete(file, allowed_exts):
    return os.path.isfile(file) and os.path.splitext(file)[1] not in allowed_exts

def main():
    parser = argparse.ArgumentParser(description="Remove all files except those with specified extensions.")
    
    parser.add_argument(
        '--keep', nargs='+', default=['.psml', '.psf', '.fdf', '.sh'],
        help="List of extensions to keep (e.g. .fdf .sh)"
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help="Only print the files that would be removed"
    )
    parser.add_argument(
        '--no-confirm', '-n', action='store_true',
        help="Do not ask for confirmation before deleting files"
    )
    parser.add_argument(
        '--path', default='.',
        help="Directory to clean (default: current directory)"
    )
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()

    allowed_exts = set(args.keep)
    
    if args.intro == True:
        show_intro()
    print("\n" + color_text("Clean:", 'bold'))
    print("-"*60)
    print("Remove all files except those with specified extensions.\n")

    for file in os.listdir(args.path):
        full_path = os.path.join(args.path, file)
        if should_delete(full_path, allowed_exts):
            if args.dry_run:
                print(f"[Dry-run] Would remove: {file}")
            elif args.no_confirm:
                os.remove(full_path)
                print(f"[INFO] Removed: {file}")
            else:
                answer = input(f"[CONFIRM] Delete {file}? [y/N] ").strip().lower()
                if answer == 'y':
                    os.remove(full_path)
                    print(f"[INFO] Removed: {file}")
                else:
                    print(f"[INFO] Skipped: {file}")
    
    print("\n[INFO] Complete job!") 
    print("\n"+"-"*60)
    print(color_text("Cleaned up! May the deleted files rest in pieces.\n\n", 'bold'))

if __name__ == "__main__":
    main()

