#!/bin/bash

set -x
set -e

python hiking.py data.txt --output data.png --figsize 6,4 --dpi 100
