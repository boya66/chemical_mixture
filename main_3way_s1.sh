#! /bin/sh 

#PBS -q batch  
#PBS -l nodes=1:R41  
#PBS -l walltime=48:00:00  
#PBS -l mem=32gb 
#PBS -m a  
#PBS -M r.lu@wustl.edu  

cd ~/project1_sensSel/

Rscript main_3way_s1.R
