#!/bin/bash


export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1




#python 5_sample_hmm_states.py BBA > bba_sample.log 2>&1 & 
#python 5_sample_hmm_states.py BBL > bbl_sample.log 2>&1 & 
#python 5_sample_hmm_states.py Chignolin > chignolin_sample.log 2>&1 & 
#python 5_sample_hmm_states.py Villin > villin_sample.log 2>&1 & 
#python 5_sample_hmm_states.py WW-domain > wwdomain_sample.log 2>&1 & 
#python 5_sample_hmm_states.py Trp-cage > trpcage_sample.log 2>&1 & 
#python 5_sample_hmm_states.py Homeodomain > homeodomain_sample.log 2>&1 & 
python 5_sample_hmm_states.py Protein-B > prb_sample.log 2>&1 & 
