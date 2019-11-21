#!/bin/bash

# number threads 
nbr_thr=( 6 )

# absolute path to input data
path_to_raw='~/fungal_ITS/'

# absolute path to output directory
path_to_data='~/fungal_ITS/'


# list the number of sequencing runs 
list=( 001 ) 

# lenght to truncate reads
len=( 0 )

# absolute path to UNITE sequence database (v7-01-12.2017) - trained classifier
class_path="~/classifier_ITS_unite_v7_011217.qza"

