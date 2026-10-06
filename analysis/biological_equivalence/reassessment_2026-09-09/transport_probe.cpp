// Focused production BioFVM probe; no FBA, growth, mechanics or original-output writes.
#include "BioFVM/BioFVM.h"
#include <omp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
using namespace BioFVM;

struct FixedAgent : Basic_Agent {
    void configure(Microenvironment& m, const std::array<double,9>& row) {
        register_microenvironment(&m);
        position = {row[0],row[1],row[2]};
        update_voxel_index();
        set_total_volume(row[3]);
        for(int d=0;d<5;++d) (*net_export_rates)[d]=row[4+d];
        set_internal_uptake_constants(.001);
    }
};

int main(int argc,char** argv) {
    // mode, input fields, input cells, output prefix, threads, steps/trials, lock, partition
    const std::string mode=argv[1], prefix=argv[4];
    const int threads=std::stoi(argv[5]), steps=std::stoi(argv[6]);
    const bool locked=std::stoi(argv[7]);
    const int partition=std::stoi(argv[8]);
    omp_set_dynamic(0); omp_set_num_threads(threads);
    int rank=0,size=1;
#ifdef PROBE_MPI
    DistPhy::mpi::mpi_Environment world{};
    DistPhy::mpi::mpi_Cartesian topo{};
    if(mode=="blocking" || mode=="avx") {
        int provided;
        MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&size);
        world.rank=rank; world.size=size; world.init_comm=MPI_COMM_WORLD;
        topo.Build_Cartesian_Topology(world);
        topo.Find_Cartesian_Coordinates(world); topo.Find_Left_Right_Neighbours(world);
    }
#endif
    Microenvironment& m=microenvironment;
    const double initial[5]={25,27.5,0,0,0};
    const double diffusion[5]={120000,35000,50000,50000,100000};
    m.set_density(0,"oxygen","dimensionless");
    m.diffusion_coefficients[0]=diffusion[0]; m.decay_rates[0]=0;
    const char* names[5]={"oxygen","glucose","acetate","lactate","CO2"};
    for(int d=1;d<5;++d) m.add_density(names[d],"dimensionless",diffusion[d],0);
#ifdef PROBE_MPI
    if(mode=="blocking" || mode=="avx")
        m.resize_space(-250.,250.,-250.,250.,-4.,4.,2.,2.,2.,world,topo);
    else
#endif
        m.resize_space(-250.,250.,-250.,250.,-4.,4.,2.,2.,2.);
    const double left=-250.+500.*rank/size, right=-250.+500.*(rank+1)/size;
#ifdef PROBE_MPI
    // Non-MPI harness modes use a full-domain local coordinate box.
    if(mode=="source" || mode=="serial") m.mesh.local_bounding_box=m.mesh.bounding_box;
#endif
    for(int n=0;n<m.number_of_voxels();++n) {
        for(int d=0;d<5;++d) m.density_vector(n)[d]=initial[d];
        const auto& pos=m.mesh.voxels[n].center;
        if(std::abs(pos[0])==249. || std::abs(pos[1])==249.) {
            std::vector<double> value(initial,initial+5);
            std::vector<bool> active(5,true);
            m.add_dirichlet_node(n,value); m.set_substrate_dirichlet_activation(n,active);
        }
    }
    default_microenvironment_options.track_internalized_substrates_in_each_agent=false;
    std::vector<FixedAgent*> agents;
    std::set<int> occupied;
    std::array<double,5> expected{};
    std::ifstream cell_input(argv[3]);
    std::array<double,9> row;
    while(cell_input>>row[0]>>row[1]>>row[2]>>row[3]>>row[4]>>row[5]>>row[6]>>row[7]>>row[8]) {
        if(row[0]<left || row[0]>=right) continue;
        if(partition==0 && row[0]>=0) continue;
        if(partition==1 && row[0]<0) continue;
        auto* a=new FixedAgent; a->configure(m,row); agents.push_back(a);
        occupied.insert(a->get_current_voxel_index());
        for(int d=0;d<5;++d) expected[d]+=.001*row[4+d];
    }
    auto sources=[&]() {
        #pragma omp parallel for
        for(int i=0;i<(int)agents.size();++i) {
            if(locked) {
                #pragma omp critical
                agents[i]->simulate_secretion_and_uptake(&m,.001);
            } else agents[i]->simulate_secretion_and_uptake(&m,.001);
        }
    };
    if(mode=="source") {
        std::ofstream out(prefix+".csv");
        out<<"trial,substrate,expected,actual,relative_error\n"<<std::setprecision(17);
        for(int trial=0;trial<steps;++trial) {
            for(int n:occupied) for(int d=0;d<5;++d) m.density_vector(n)[d]=initial[d];
            sources();
            for(int d=0;d<5;++d) {
                double actual=0;
                for(int n:occupied) actual+=(m.density_vector(n)[d]-initial[d])*m.mesh.voxels[n].volume;
                out<<trial<<','<<names[d]<<','<<expected[d]<<','<<actual<<','
                   <<(expected[d] ? (actual-expected[d])/expected[d] : 0)<<'\n';
            }
        }
        std::cout<<"agents="<<agents.size()<<" occupied_voxels="<<occupied.size()
                 <<" threads="<<threads<<" locked="<<locked<<" trials="<<steps<<'\n';
        return 0;
    }
    std::ifstream input(argv[2],std::ios::binary);
    std::array<double,8> data;
    while(input.read(reinterpret_cast<char*>(data.data()),sizeof(data))) {
        if(data[0]<left || data[0]>=right) continue;
        std::vector<double> pos={data[0],data[1],data[2]};
        int n=m.nearest_voxel_index(pos);
        for(int d=0;d<5;++d) m.density_vector(n)[d]=data[3+d];
    }
    for(int step=0;step<steps;++step) {
#ifdef PROBE_MPI
        if(mode=="blocking") diffusion_decay_solver__constant_coefficients_LOD_3D_BLOCKING(m,.001,world,topo);
        else if(mode=="avx") diffusion_decay_solver__constant_coefficients_LOD_3D_AVX256D(m,.001,world,topo);
        else
#endif
            diffusion_decay_solver__constant_coefficients_LOD_3D(m,.001);
        if(locked) sources(); // Fixed-source replay with one thread; no FBA calls.
    }
    std::ofstream output(prefix+"_"+std::to_string(rank)+".bin",std::ios::binary);
    for(int n=0;n<m.number_of_voxels();++n) {
        for(int d=0;d<3;++d) data[d]=m.mesh.voxels[n].center[d];
        for(int d=0;d<5;++d) data[d+3]=m.density_vector(n)[d];
        output.write(reinterpret_cast<char*>(data.data()),sizeof(data));
    }
    std::cout<<"mode="<<mode<<" rank="<<rank<<" steps="<<steps<<" source_replay="<<locked<<'\n';
#ifdef PROBE_MPI
    if(mode=="blocking" || mode=="avx") MPI_Finalize();
#endif
}
