# Iliad - General Relativistic Pulsar Timing in Kerr Spacetime

Welcome to Iliad! This code can be used to model the time-frequency signal from a pulsar orbiting a massive black hole (BH), accounting for all relativistic and astrophysical effects. Some work using these methods has been published in [Kimpson 2019a](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2411K/abstract),[b](https://ui.adsabs.harvard.edu/abs/2019MNRAS.486..360K/abstract). 

**The public release of this code via this repo is still in development.** Whilst the bare-bones base code is available for constructing the signal, modules relating to the signal analysis and computational optimisation are still due to be released, once they have been cleaned up!


This repo  has two main ingredients. The first part determines the orbital trajectory of a spinning body in a curved spacetime. It is based on the [Spin Curvature Dynamics](https://github.com/tomkimpson/SpinCurvatureDynamics) code. The second part calculates the trajectory of light from the pulsar, based on the [ForwardRayTracing](https://github.com/tomkimpson/ForwardRayTracing) code. This repo then combines these two tools so that the user can specify some BH-PSR system and consistently generate the frequency-dependent photon ToAs. 



## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

The main part of the code (w.r.t orbital dynamics, ray tracing) is written in FORTRAN with a [gfortran](https://gcc.gnu.org/wiki/GFortran) compiler. **Other compilers have not been tested.** The gfortran installation binaries can be found [here](https://gcc.gnu.org/wiki/GFortranBinariels), although typically gfortran comes pre-installed on most Linux/Unix systems. If you have [Homebrew](https://brew.sh/) installed on OSX, you can simply run 


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

You can check the environment variable has been added to `bash_profile` by either `env` or `vim ~/.bashprofile`


The code should then run as is, out of the box. Try

```
run.py
```
to compile and run the code. Once you have checked that everything is running OK, you can then start playing.The `Worked_Example` Jupyter notebook is a good place to start, for a python interface with the underlying code. The code structure (modules, subroutines etc.) is outlined below.


If making edits to the code, try to keep to the [FORTRAN Style Guide](https://www.fortran90.org/src/best-practices.html)


## Structure

* `parameters.f`. This is effectively the config file for specifying the system parameters

* `constants.f`. This is a globally accessible module that performs some calculations at compile time based on the specifications in `parameters.f`. Nothing should need to be altered in this file.

* `main.f`. This is the root program which runs all the necessary subroutines. It first determines the PSR spin orbital dynamics (`OrbitalDynamics.f`) and then uses the output as initial conditions for the ray tracing (`RayTracing.f`) 

* `OrbitalDynamics.f`. Computes the PSR spin dynamics, assuming a Kerr background spacetime. The initial position, momenta and spin can all be set in `parameters.f`. Similarly, the numerical integration method can be set to have a constant or adaptive stepsize.The output is a single binary file of size NSTEPS x 13 i.e. number of integration steps and 3 x 4 vectors (positions, spin ,momentum) plus the proper time (related to the stepsize). If `plot_MPD` is set to 1 in `parameters.f`, then the output file will also be written as a readable  text file (useful for e.g. plotting).

* `RayTracing.f`. Computes the path of the pulsar light ray on the background spacetime. Either a 'forward' or 'backward' in time integration can be used, depending on your goals (e.g. TOAs vs examination of precession, pulse profile etc.). It takes the output of `OrbitalDynamics` as boundary conditions and finds the correct ray using a non-linear conjugate gradient descent approach.

* `NumericalMethods.f`. Governs the numerical integration. We use a 5th order Runge-Kutta with adaptive stepsize.

* `derivatives`. Defines the derivative ODEs to be solved numerically.

* `metric.f`. Completely defines the background spacetime and associated Christoffel, Riemann etc. components.

* `IO.f`. Module for IO. Typically unformatted binary. Exploring use of HDF.



![](example_data/flowchart.png?raw=true)




##


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details




