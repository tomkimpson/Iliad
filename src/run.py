from os import system as os


os("gfortran -fopenmp -ffree-form -fdefault-real-8 -O3 *.f -J mod/") 
os("gfortran -fopenmp -ffree-form -fdefault-real-8 -O3 *.f -J mod/") 

os("./a.out")
