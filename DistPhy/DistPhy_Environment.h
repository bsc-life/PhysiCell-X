#include<mpi.h>

#ifndef __DistPhy_Environment_h__
#define __DistPhy_Environment_h__

namespace DistPhy
{
    namespace mpi
    {
        class mpi_Environment
        {
            public:
                int         size;               //Size of MPI_COMM_WORLD
                int         rank;               //Rank of each MPI process
                MPI_Comm    init_comm;          //Initial communicator i.e. MPI_COMM_WORLD
              
            public:
                void Initialize();              //Results in MPI_Init_thread(NULL, NULL, req_thread_level=MPI_THREAD_MULTIPLE,&provided_thread_level);
                int  Size();                    //Returns size of Communicator "comm"
                int  Rank();                    //Returns rank of MPI process in communicator "comm"
                void Finalize();                //Simply calls MPI_Finalize();             
        };
        
        //Non-class functions in namespace DistPhy::mpi
        int IOProcessor(mpi_Environment &world);
       
    }
}

#endif
        
    
