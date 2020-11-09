DistPhy (Distributed PhysiCell) is an object-oriented parallel framework for use in PhysiCell/BioFVM. All entities within DistPhy reside within the namespace "DistPhy". 
A second level namespace "mpi" encapsulates MPI specific functionalities. 
The following is the aim:
(a) DistPhy must not be dependent on other libraries except for MPI.
(b) Interfaces in DistPhy must be as simple as possible. 
(c) Functions of DistPhy should abstract as many details as possible from the user.
(d) Wherever possible DistPhy must maintain an object-oriented spirit.