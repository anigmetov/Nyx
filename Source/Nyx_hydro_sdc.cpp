
#include "Nyx.H"
#include "Nyx_F.H"
#include <AMReX_Particles_F.H>

#ifdef GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;

using std::string;

#ifndef NO_HYDRO
void
Nyx::sdc_hydro (Real time,
                Real dt,
                Real a_old,
                Real a_new)
{
    BL_PROFILE("Nyx::sdc_hydro()");

    BL_ASSERT(NUM_GROW == 4);

    int sdc_iter;
    const Real prev_time    = state[State_Type].prevTime();
    const Real cur_time     = state[State_Type].curTime();

    MultiFab&  S_old        = get_old_data(State_Type);
    MultiFab&  S_new        = get_new_data(State_Type);

    MultiFab&  D_old        = get_old_data(DiagEOS_Type);
    MultiFab&  D_new        = get_new_data(DiagEOS_Type);

    MultiFab&  IR_old       = get_old_data(SDC_IR_Type);
    MultiFab&  IR_new       = get_new_data(SDC_IR_Type);

    // This is just a hack to get us started
    if (time == 0.)
    {
       IR_old.setVal(0.);
       IR_new.setVal(0.);
    }

    // It's possible for interpolation to create very small negative values for
    // species so we make sure here that all species are non-negative after this
    // point
    enforce_nonnegative_species(S_old);

    MultiFab ext_src_old(grids, dmap, NUM_STATE, 3);
    ext_src_old.setVal(0.);

    // Define the gravity vector
    MultiFab grav_vector(grids, dmap, BL_SPACEDIM, 3);
    grav_vector.setVal(0.);

#ifdef GRAVITY
    gravity->get_old_grav_vector(level, grav_vector, time);
    grav_vector.FillBoundary(geom.periodicity());
#endif

    // Create FAB for extended grid values (including boundaries) and fill.
    MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
    FillPatch(*this, S_old_tmp, NUM_GROW, time, State_Type, 0, NUM_STATE);

    MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);
    FillPatch(*this, D_old_tmp, NUM_GROW, time, DiagEOS_Type, 0, D_old.nComp());

    MultiFab hydro_src(grids, dmap, NUM_STATE, 0);
    hydro_src.setVal(0.);

    MultiFab divu_cc(grids, dmap, 1, 0);
    divu_cc.setVal(0.);

    //Begin loop over SDC iterations
    int sdc_iter_max = 2;

    for (sdc_iter = 0; sdc_iter < sdc_iter_max; sdc_iter++)
    {
       amrex::Print() << "STARTING SDC_ITER LOOP " << sdc_iter << std::endl;

       // FillPatch the IR term into ext_src_old for both Eint and Eden
       // ASA -- Is this term actually the source term for (rho e) or (e)???
       MultiFab IR_tmp(grids, dmap, 1, NUM_GROW);
       FillPatch(*this, IR_tmp, NUM_GROW, time, SDC_IR_Type, 0, 1);
       MultiFab::Copy(ext_src_old,IR_tmp,0,Eden,1,0);
       MultiFab::Copy(ext_src_old,IR_tmp,0,Eint,1,0);

       bool   init_flux_register = (sdc_iter == 0);
       bool add_to_flux_register = (sdc_iter == sdc_iter_max-1);
       compute_hydro_sources(time,dt,a_old,a_new,S_old_tmp,D_old_tmp,
                             ext_src_old,hydro_src,grav_vector,divu_cc,
                             init_flux_register, add_to_flux_register);

       // NOTE: WE DON'T ACTUALLY WANT TO UPDATE THE WHOLE STATE HERE
       //       THIS IS JUST WHAT WE DO IN THE STRANG VERSION
       update_state_with_sources(S_old_tmp,S_new,
                                 ext_src_old,hydro_src,grav_vector,divu_cc,
                                 dt,a_old,a_new);

       // We copy old Temp and Ne to new Temp and Ne so that they can be used
       //    as guesses when we next need them.
       MultiFab::Copy(D_new,D_old,0,0,D_old.nComp(),0);
       
       // ASA -- Why do we do this??  Did we change S_old_tmp at some point
       MultiFab::Copy(S_old_tmp,S_old,0,0,S_new.nComp(),0);

       // This step needs to do the update of (rho),  (rho e) and (rho E)
       //      AND  needs to return an updated value of I_R in the old SDC_IR statedata.
       sdc_reactions(S_old_tmp, S_new, D_new, hydro_src, IR_old, dt, a_old, a_new, sdc_iter);

       // ASA -- This is just a hack for now to see what happens
       IR_old.setVal(0.);

    } //End loop over SDC iterations

    // Copy IR_old (the current IR) into IR_new here so that when the pointer swap occurs
    //     we can use IR in the next timestep
    MultiFab::Copy(IR_new,IR_old,0,0,1,0);

    grav_vector.clear();
}
#endif
