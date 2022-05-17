#!/bin/bash


export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1


#python 6_CK_test.py Homeodomain > homeodomain_cktest.log 2>&1 & 
python 6_CK_test.py WW-domain > wwdomain_cktest.log 2>&1 & 
#python 6_CK_test.py BBL > bbl_cktest.log 2>&1 & 
#python 6_CK_test.py Protein-B > prb_cktest.log 2>&1 & 
