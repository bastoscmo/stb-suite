#!/bin/bash
for d in strain_*/ ; do
    [ -d "$d" ] && echo "Checking $d..." && tail -n 1 "$d/calc.out" 2>/dev/null
done
