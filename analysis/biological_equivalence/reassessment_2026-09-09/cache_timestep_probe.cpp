// Reproduce the dFBA-to-BioFVM cached timestep handoff with production methods.
#include "BioFVM/BioFVM.h"
#include <iostream>
#include <iomanip>
using namespace BioFVM;
int main() {
    auto& m=microenvironment;
    m.resize_space(-1.,1.,-1.,1.,-1.,1.,2.,2.,2.);
#ifdef PROBE_MPI
    m.mesh.local_bounding_box=m.mesh.bounding_box;
#endif
    default_microenvironment_options.track_internalized_substrates_in_each_agent=false;
    Basic_Agent a;
    a.register_microenvironment(&m);
    a.position={0.,0.,0.}; a.update_voxel_index();
    a.set_total_volume(1.);
    (*a.net_export_rates)[0]=8.;
    m.density_vector(0)[0]=0.;
    // update_dfba_outputs builds cached increments using intracellular_dt=0.01.
    a.set_internal_uptake_constants(.01);
    // Secretion::advance calls this every diffusion_dt=0.001 without a cache rebuild.
    a.simulate_secretion_and_uptake(&m,.001);
    std::cout<<std::setprecision(17)<<"cached_dt=0.01 call_dt=0.001 expected=0.001 actual="
             <<m.density_vector(0)[0]<<'\n';
    m.density_vector(0)[0]=0.;
    a.set_internal_uptake_constants(.001);
    a.simulate_secretion_and_uptake(&m,.001);
    std::cout<<"cached_dt=0.001 call_dt=0.001 expected=0.001 actual="
             <<m.density_vector(0)[0]<<'\n';
}
