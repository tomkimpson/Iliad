from os import system as os

#Clear data directories in preparation for fresh run
#os("rm $IliadDir*.dat")
#os("rm $IliadDir*.txt")


#also clear the example_data directories

#os("rm ../example_data/*.dat")
#os("rm ../example_data/*.txt")



#Compile the code. 
#Compilation flags include -fopenmp for OpenMP parallelisation and -Werror to prevent compilation if warning messages are present
os("gfortran -fopenmp -Werror -ffree-form -fdefault-real-8 -Ofast parameters.f constants.f IO.f metric.f derivatives.f NumericalMethods.f OrbitalDynamics.f RayTracing.f main.f -J mod/") 

#run it
os("./a.out")
