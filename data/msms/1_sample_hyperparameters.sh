#!/bin/zsh

#  -h, --help            show this help message and exit
#  -i SEARCH_SPACE, --search-space SEARCH_SPACE
#                        Path to file that defines the HP search space dictionary
#  -n NUM_TRIALS, --num-trials NUM_TRIALS
#                        Number of HP trials
#  -o OUTPUT_FILE, --output-file OUTPUT_FILE
#                        Path to hd5 file to store HP samples
#  -s SEED, --seed SEED  Random seed

#conda activate msmsense

msmsense sample \
 -n 10000 \
 -o hpsample_long.h5


# RUN THE TRUNCATE NOTEBOOK TO GET AN EVEN SAMPLE ACROSS ALL DIMENSIONS.