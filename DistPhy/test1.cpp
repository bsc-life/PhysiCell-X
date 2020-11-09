#include<iostream>

#include"DistPhy.h"

using namespace std;
using namespace DistPhy::mpi; 

int main()
{
    int my_rank;
    
    mpi_Environment world; 
    world.Initialize(); 
    
    mpi_Cartesian cart; 
    cart.Build_Cartesian_Topology(world);
    
    if(IOProcessor(world))
        cout<<"World Size = " << world.size; 
    
    my_rank = world.Rank(); 
    
    world.Finalize(); 
}
