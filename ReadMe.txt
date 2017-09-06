This code is a path integral Monte Carlo code developed by
Yangqian Yan. The majority of the code was developed during the
Ph.D. studies of the author, conducted in the group of
Prof. Doerte Blume at Washington State University. 


Use of the code is at the risk of the user.

The code is an out-growth of the codes that led to the following
publication:

Phys. Rev. Lett. 116, 230401 (2016), "Path-Integral Monte Carlo
Determination of the Fourth-Order Virial Coefficient for a Unitary
Two-Component Fermi Gas with Zero-Range Interactions" by
Yangqian Yan and D. Blume.

Related technical details can be found in the following
publications:

Ref.1 Manuscript currently under review with the Journal of Physics B,
"Path integral Monte Carlo ground state
approach: Formalism, implementation, and applications" by
Yangqian Yan and D. Blume.

Ref.2 Phys. Rev. A 91, 043607 (2015), "Incorporating exact two-body
propagators for zero-range interactions into N -body Monte Carlo
simulations" by Yangqian Yan and D. Blume.

and

Ref.3 Ph.D. Thesis entitled "Path integral Monte Carlo studies of
ultracold few-atom systems" by Yangqian Yan, Washington State
University (Washington State University, ProQuest Dissertations
Publishing, 2016. 10139760).

Support by the NSF, which made this work possible, is gratefully
acknowledged.


The code is a mixture of C and C++. 
The scalar estimators and wiggle/permute moves are written in C.
The whole path and center of mass and dilation moves, as well as the correlation
estimators, are written in C++.

==============================
Requirements:
==============================
0. Linux environment with shell access.
1. C/C++ compiler: including GNU or INTEL
2. The Message Passing Interface (MPI): MPICH or Open MPI
3. GSL - GNU Scientific Library

Depending on the complier and MPI configuration, the wrapper might be named
differently.
If Open MPI and the gnu compiler are used, the wrapper is mpic++.
The code is tested on
C++03: g++ (Ubuntu 4.8.4-2ubuntu1~14.04) 4.8.4 + mpirun (Open MPI) 1.10.2
C++11: g++ (GCC) 6.1.1 20160501 + mpirun (Open MPI) 1.10.2



==============================
Compile and execute
==============================
To compile the code, run "make".
Please make sure the gsl library is pointed to correctly. 
The resulting objects and executable will be in the directory "fermi".

Run "make clean" to clean all objects, executables, and test results.

Run "run --help" to get a list of the options that the program takes.
See "job.sh" for a bash script to run the program.
To run the code on a computer cluster, the script needs to be modified 
accordingly.

==============================
Output explanation
==============================
The results and various comments are directly written to the screen.
It is recommended that this output be piped into a file.
For scalar estimators, the output consists of the name of the estimator,
the estimated mean, the estimated error, and the statistical error without the 
sign error.

The data can be further processed using tools such as grep, sed, awk, or
using python/perl (or any other tool).
Detailed data for the first scalar estimator is recorded in "col0.dat" during
the simulation. "col0.dat" is formatted such that each row contains one step of
the simulation. The first column is the estimate of the sign and the remaining
columns contain the observables.
One can monitor the result by running "cat col0.dat|b4” while the simulation 
is still in progress.
Here, “b4” is a bash script located in the root directory.
The resulting file contains 80*N measurements, where N is the number of 
processors.
The output consists of the average, the error due to the statistics, 
the error due to the sign, and the total error.
This is the preferred way of measuring the error.

The error that is printed directly by the program to the screen (treating N processors
as N measurements) is calculated differently than the bash script “b4”. 
For a large number of processors, the two should roughly agree. For a small
number of processors, the error printed directly by the program to the 
screen is not as accurate.

The data for the second observable are stored in "col1.dat", etc.

The scaled structural properties are written to the file "fwave#.dat".
The unscaled structural properties are written to the file "fwave#d.dat".
E.g.: the scaled pair distribution is stored as "fwave2.dat" and
the unscaled single particle distribution is stored as "fwave1d.dat".
Typical header reads "x average std 0 1 2", which means 
position, histogram average, standard deviation, distribution 
for particle 0, particle 1, and particle 2.
The mean (column 2) and standard deviation (column 3) are calculated using 
columns 4,5, and 6.


==============================
Tests
==============================
The directory "./test/" contains an example script (named "job.sh")
for a simulation of 2 particles and the
corresponding sample output (file name "sample.out").
The file "job.sh" contains comments and should be selfexplanatory.
The output should be identical to "sample.out" except for the time and date.

The directory "./test/3” contains another example script (named "job.sh") and
the corresponding sample output (file named "sample.out").
This example is for 3 particles and takes longer to run. The random number seed
is also not fixed. This example can be fairly easily changed to an actual production
run/simulation.

