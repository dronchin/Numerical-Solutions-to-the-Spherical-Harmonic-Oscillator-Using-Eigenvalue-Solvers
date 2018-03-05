# Project2
PHY 480 with professor Morten Hjorth

Program source files can be found in the program source folder. This includes catch.hpp, EigenvalueArma.cpp, functiontest.cpp, jacobiMethod.cpp, jacobiMethod.h, justjacobi.cpp, main.cpp, mainMyjacobi.cpp, plottingDifEnergy.py, plottinginteraction.py

Use the function “make” in the Source files folder in order to build the three programs.
Builds main.x, functiontest.x, and mevsArma.x files. All programs were compiled using g++ and the -O3  flag as well as the following libraries: -larmadillo -llapack -lblas.

For the main.x, it will ask you to input N which will determine the size of the matrix as nxn. Pick an n between 100-1000. The larger values of N take longer to run but provide more numerical accuracy. Next, it will ask for a omega value. This determines the strength of the harmonic oscillator potential. It will then ask for a rho_max. This will determine what energy values are calculated and the precision of those numbers. Pick a rho_max so that u(rho_max) ~= 0. Finally, it will ask if you want the interacting or non-interacting case. “Y” or “y” for the interacting case and anything else for non-interacting case. This program will write 3 files and print out the first 5 eigenvalues. The 3 files are needed for the python plotting files.

For functiontest.x, no inputs are required. All tests are run using Catch.cpp. This project can be found at https://github.com/catchorg/Catch2

For mevsArma.x, It will ask for the same inputs as main.x, I would recommend a lower value for N in the 50-400 range as it will run my Jacobi Eigenvalue solver which takes more time compared to the armadillo eig_sym() funciton.

The report source code and pdf can be found in the Report file.
