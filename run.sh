#!/bin/bash

set -x
set -e

python hiking.py data.txt --output data.png --figsize 8,6 --dpi 100
