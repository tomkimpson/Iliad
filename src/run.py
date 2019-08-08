from os import system as os

#Clear data directory

os("rm $IliadDir/*.dat")

os("gfortran -fopenmp -Werror -ffree-form -fdefault-real-8 -O3 *.f -J mod/") 
os("gfortran -fopenmp -Werror -ffree-form -fdefault-real-8 -O3 *.f -J mod/") 

os("./a.out")
