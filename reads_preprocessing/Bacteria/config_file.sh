#!/bin/bash

# number threads 
nbr_thr=( 6 )

# absolute path to input data
path_to_raw='~/MWDW'

# absolute path to output directory
path_to_data='~/MWDW'


# list the number of sequencing runs 
list=( 001 ) 

# lenght to truncate reads
len=( 0 )

# absolute path to Green Genes sequence database (v08.13) - trained classifier on 997F-1192R
class_path="~/gg-13-8-99-977F-1192R-classifier.qza"

