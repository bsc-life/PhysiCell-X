# PhysiCell-X: A Distributed-Shared Parallel Version of PhysiCell

**Version:** 1.14

**Release date:** 15 March 2026

**Reference:**

	(a) The PhysiCell-X repository
	
	(b) Saxena, Gaurav, Miguel Ponce-de-Leon, Arnau Montagud, David Vicente Dorca, and Alfonso Valencia. 
	"BioFVM-X: An MPI+ OpenMP 3-D Simulator for Biological Systems."
	In International Conference on Computational Methods in Systems Biology, pp. 266-279. Springer, Cham, 2021.

	(c) Jose-Luis Estragues-Munoz, Carlos Alvarez, Arnau Montagud, Daniel Jimenez-Gonzalez, and Alfonso Valencia. 
	"A Novel Scalable High Performance diffusion solver for multiscale cell simulations."
	In preprint

**User Guide:** PhysiCell-X_UserGuide.pdf in the documentation folder (use in conjunction with the PhysiCell documentation). 


PhysiCell-X is the distributed-shared parallel version of PhysiCell (http://physicell.org). PhysiCell is a open-source,
multi-physics, multi-scale, agent-based simulator for biological systems. It provides both the stage
(micro-environment) and actors (cells or agents) for simulation. Though PhysiCell is light-weight, flexible and
shared-memory parallelized using OpenMP (Open Multiprocessing), it cannot run on distributed systems i.e.
it cannot be executed on multiple nodes of an HPC (High Performance Computing) cluster. Thus, the problem size that
PhysiCell can handle is limited by the maximum memory of a single node. This is a limitation that needs to be
removed and this is where PhysiCell-X comes in. PhysiCell-X enables the distributed parallelization of PhysiCell by
making use of MPI (Message-Passing Interface). In simple words, you can now use multiple nodes of an HPC system to
solve a single, coherent problem. Thus, the aim of PhysiCell-X is to remove the memory limitation, reduce the time
to solution by distributing the computation onto multiple compute nodes and solve very large sized (real-world
scale) problem.

Introduction
============

PhysiCell-X is based on PhysiCell v1.9.0 which also added support for intracellular modelling. The current version of 
PhysiCell-X is v0.1 and it is actively being developed. PhysiCell-X as of now supports only 3-D problems.

Pre-requisites
==============
There are two pre-requisites to running/modelling projects with PhysiCell-X:

	(a) A C++ compiler for e.g. GNU compiler (preferred), Intel compiler etc. that supports OpenMP
	(b) An MPI implementation such as OpenMPI (preferred), Intel MPI, MPICH.  

Sample Projects
===============
There are multiple MPI parallelized sample projects that can be found in the directories **sample_projects** 
and **sample_projects_intracelluar** . To see a list of available projects, do  

`$ make list-projects`

The projects which are MPI parallelized (PhysiCell-X) have *-mpi* as the suffix. To check out a certain project (like
heterogeneity-sample-mpi), do

`$ make heterogeneity-sample-mpi` 

copies multiple files in generic directories under the top-level directory like **./config, ./custom_modules.** 
To compile this checked out project, do

`$ make` 

After the executable is produced, it can be submitted for execution using one of the provided submission script files (specifically
for SLURM) or executed simply using `mpiexec` or `mpirun`. 

The working directories can be cleaned using 

`$ make clean` 

and reset to the original state by executing

`$ make reset`

There are other Make rules which can be seen in the top level `Makefile`. 

PhysiCell (http://physicell.org)
================================
Pure OpenMP programs inside PhysiCell-X can also be executed as PhysiCell-X v0.1 is based on PhysiCell-v1.9.0 but
one needs to compile the pure OpenMP examples through `mpic++` i.e. the MPI compiler wrapper.

Important
=========
Please see the detailed user documentation titled **PhysiCell-X_UserGuide.pdf** inside the documentation directory. 



