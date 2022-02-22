#!/bin/bash

python model_check.py 1fme 0 > 0.log 2>&1 &
python model_check.py 1fme 1 > 1.log 2>&1 &
python model_check.py 1fme 4 > 4.log 2>&1 &