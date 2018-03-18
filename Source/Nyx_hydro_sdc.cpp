
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

    int sdc_iter;
    const Real prev_time    = state[State_Type].prevTime();
    const Real cur_time     = state[State_Type].curTime();
    MultiFab&  S_old        = get_old_data(State_Type);
    MultiFab&  D_old        = get_old_data(DiagEOS_Type);
    MultiFab&  S_new        = get_new_data(State_Type);
    MultiFab&  D_new        = get_new_data(DiagEOS_Type);

    if (std::abs(time-prev_time) > (1.e-10*cur_time) )
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "sdc_hydro:  prev_time = " << prev_time << std::endl;
            std::cout << "sdc_hydro:       time = " <<      time << std::endl;
        }
        amrex::Abort("time should equal prev_time in sdc_hydro!");
    }

    // It's possible for interpolation to create very small negative values for
    // species so we make sure here that all species are non-negative after this
    // point
    enforce_nonnegative_species(S_old);

    MultiFab ext_src_old(grids, dmap, NUM_STATE, 3);
    ext_src_old.setVal(0);

    //Use previous I_R, ignore growth cells
    // TREATING I_R as zero
    MultiFab I_R(grids, dmap, 1, 3);
    I_R.setVal(0.);
    // MultiFab::Copy(ext_src_old,D_old,Diag1_comp,Eint,1,0);
    
    // Define the gravity vector so we can pass this to ca_umdrv.
    MultiFab grav_vector(grids, dmap, BL_SPACEDIM, 3);
    grav_vector.setVal(0.);

#ifdef GRAVITY
    gravity->get_old_grav_vector(level, grav_vector, time);
    grav_vector.FillBoundary(geom.periodicity());
#endif

    BL_ASSERT(NUM_GROW == 4);

    // Create FAB for extended grid values (including boundaries) and fill.
    MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
    FillPatch(*this, S_old_tmp, NUM_GROW, time, State_Type, 0, NUM_STATE);

    MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);
    FillPatch(*this, D_old_tmp, NUM_GROW, time, DiagEOS_Type, 0, D_old.nComp());

    MultiFab hydro_src(grids, dmap, NUM_STATE, 0);
    hydro_src.setVal(0);

    MultiFab divu_cc(grids, dmap, 1, 0);
    divu_cc.setVal(0);

    FArrayBox flux[BL_SPACEDIM], u_gdnv[BL_SPACEDIM];

    //Begin loop over SDC iterations
    int sdc_iter_max = 2;

    for (sdc_iter = 0; sdc_iter < sdc_iter_max; sdc_iter++)
    {
       amrex::Print() << "STARTING SDC_ITER LOOP " << sdc_iter << std::endl;

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
       
       MultiFab::Copy(S_old_tmp,S_old,0,0,S_new.nComp(),0);

         // Gives us I^1 from integration with F^(n+1/2) source
         // Stores I^1 in D_new(diag1_comp) and ext_src_old(UEINT)

       sdc_reaction_step(time, dt, S_old_tmp, D_new, ext_src_old,sdc_iter);

	 // Consider changing to a copy operation
	 // I_R is stored in two places, but we need it in diag_eos for the next timestep
	 //	 MultiFab::Copy(ext_src_old,D_new,Diag1_comp,Eint,1,0);
    
	 //If another sdc iteration is performed, use the old state data
	 // This could be more general (sdc_iter< sdc_iter_max-1)
        if (sdc_iter < sdc_iter_max-1) 
        {
	    ext_src_old.setVal(0);
	    MultiFab::Copy(ext_src_old,I_R,0,Eint,1,0);
	    MultiFab::Copy(S_old_tmp,S_old,0,0,S_old.nComp(),0);

        } else {

	  printf("Density component is: %i",Density);
	  MultiFab::Copy(S_new,S_old_tmp,Density,Density,1,0);

	  /////////////When adding src for rho_e and rho_E works in integrate state, uncomment 
	  // these lines:
	  //  	  MultiFab::Copy(S_new,S_old_tmp,Eden,Eden,1,0);
	  //	  MultiFab::Copy(S_new,S_old_tmp,Eint,Eint,1,0);
	  //MultiFab::Copy(D_new,D_old,Temp_comp,Temp_comp,1,0);
	}

    }
    //End loop over SDC iterations

    grav_vector.clear();

}
#endif
