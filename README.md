# Iliad - General Relativistic Pulsar Timing in Kerr Spacetime

This code calculates the time-frequency behaviour from a pulsar orbiting a massive black hole (BH).

This code has two main ingredients. The first part determines the orbital trajectory of a spinning body. It is based on the [Spin Curvature Dynamics](https://github.com/tomkimpson/SpinCurvatureDynamics) code. The second part calculates the trajectory of light from the pulsar, based on the [ForwardRayTracing](https://github.com/tomkimpson/ForwardRayTracing) code. This repo combines these two tools so that the user can specify some BH-PSR system and consistently generate the frequency-dependent photon ToAs.


## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

This code is written in FORTRAN with a [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. **Other compilers have not been tested.** The gfortran installation binaries can be found [here](https://gcc.gnu.org/wiki/GFortranBinariels), although typically gfortran comes pre-installed on most Linux/Unix systems. If you have [Homebew](https://brew.sh/) installed on OSX, you can simply run 


```
brew install gcc
```

### First steps
After [cloning the repo](https://help.github.com/en/articles/cloning-a-repository), the first thing to do is to set the path to the output files that the code will produce.
This can be done by setting the environment variable as

```
echo 'export IliadDir="/Users/tomkimpson/Data/Iliad/"' >> ~/.bash_profile
source ~/.bash_profile
```

Just change the path `Users/tomkimpson/Data/Iliad/` to some appropriate local path. 

You can check the environemnt variable has been added to `bash_profile` by either `env` or `vim ~/.bashprofile`


The code should then run as is, out of the box. Try

```
run.py
```
to compile and run the code. Once you have checked that everything is running OK, you can then start playing. The code structure (modukes, subroutines etc.) is outlined below.


If making edits to the code, try to keep to the [FORTRAN Style Guide](https://www.fortran90.org/src/best-practices.html)


## Structure

* `main.f`. This is the root program which runs all the necessary subroutines. It first determines the PSR spin orbital dynamics (`OrbitalDynamics.f`) and then uses the output as initial conditions for the ray tracing (`RayTracing.f`) 

##

