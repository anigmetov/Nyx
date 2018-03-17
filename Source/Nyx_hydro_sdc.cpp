
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
    const int  finest_level = parent->finestLevel();
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

    if (do_reflux && level < finest_level)
    {
        //
        // Set reflux registers to zero.
        //
        get_flux_reg(level+1).setVal(0);
    }
    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister* fine    = 0;
    FluxRegister* current = 0;

    if (do_reflux && level < finest_level)
        fine = &get_flux_reg(level+1);
    if (do_reflux && level > 0)
        current = &get_flux_reg(level);

    const Real* dx     = geom.CellSize();
    Real        courno = -1.0e+200;

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

    MultiFab fluxes[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        fluxes[j].define(getEdgeBoxArray(j), dmap, NUM_STATE, 0);
        fluxes[j].setVal(0.0);
    }

    BL_ASSERT(NUM_GROW == 4);

    Real  e_added = 0;
    Real ke_added = 0;

    // Create FAB for extended grid values (including boundaries) and fill.
    MultiFab S_old_tmp(S_old.boxArray(), S_old.DistributionMap(), NUM_STATE, NUM_GROW);
    FillPatch(*this, S_old_tmp, NUM_GROW, time, State_Type, 0, NUM_STATE);

    MultiFab D_old_tmp(D_old.boxArray(), D_old.DistributionMap(), D_old.nComp(), NUM_GROW);
    FillPatch(*this, D_old_tmp, NUM_GROW, time, DiagEOS_Type, 0, D_old.nComp());

    FArrayBox flux[BL_SPACEDIM], u_gdnv[BL_SPACEDIM];

    //Begin loop over SDC iterations
    int sdc_iter_max = 2;

    for (sdc_iter = 0; sdc_iter < sdc_iter_max; sdc_iter++)
    {
       amrex::Print() << "STARTING SDC_ITER LOOP " << sdc_iter << std::endl;
      
#ifdef _OPENMP
#pragma omp parallel reduction(max:courno) reduction(+:e_added,ke_added)
#endif
          {
          Real cflLoc = -1.e+200;

          for (MFIter mfi(S_old_tmp,true); mfi.isValid(); ++mfi)
          {
	    const Box& bx        = mfi.tilebox();

             // Allocate fabs for fluxes.
             for (int i = 0; i < BL_SPACEDIM ; i++) {
                 const Box &bxtmp = amrex::surroundingNodes(bx, i);
                 flux[i].resize(bxtmp, NUM_STATE);
                 u_gdnv[i].resize(amrex::grow(bxtmp, 1), 1);
                 u_gdnv[i].setVal(1.e200);
             }

             FArrayBox& S_state    = S_old_tmp[mfi];
             FArrayBox& S_stateout = S_new[mfi];

#ifdef SHEAR_IMPROVED
             FArrayBox& am_tmp = AveMom_tmp[mfi];
#endif

             Real se  = 0;
             Real ske = 0;

  	    // Get F^(n+1/2)

	fort_advance_gas
            (&time, bx.loVect(), bx.hiVect(), 
             BL_TO_FORTRAN(S_state),
             BL_TO_FORTRAN(S_stateout),
             BL_TO_FORTRAN(u_gdnv[0]),
             BL_TO_FORTRAN(u_gdnv[1]),
             BL_TO_FORTRAN(u_gdnv[2]),
             BL_TO_FORTRAN(ext_src_old[mfi]),
             BL_TO_FORTRAN(grav_vector[mfi]),
             dx, &dt,
             BL_TO_FORTRAN(flux[0]),
             BL_TO_FORTRAN(flux[1]),
             BL_TO_FORTRAN(flux[2]),
             &cflLoc, &a_old, &a_new, &se, &ske, 
             &print_fortran_warnings, &do_grav, &sdc_split);
	
        for (int i = 0; i < BL_SPACEDIM; ++i) 
          fluxes[i][mfi].copy(flux[i], mfi.nodaltilebox(i));

         e_added += se;
        ke_added += ske;
       } // end of MFIter loop

        courno = std::max(courno, cflLoc);

       } // end of omp parallel region

       // If at last iteration of sdc, now have fluxes in stateout (S_new) as well as momentum etc
       
       // We copy old Temp and Ne to new Temp and Ne so that they can be used
       //    as guesses when we next need them.
       MultiFab::Copy(D_new,D_old,0,0,D_old.nComp(),0);
       
       //MultiFab::Copy(S_old_tmp,S_new,Eden,Eden,1,0);
       /*
       MultiFab::Copy(S_old_tmp,S_new,Eint,Eint,1,0);
       MultiFab::Copy(S_old_tmp,S_new,0,0,S_new.nComp(),0);*/
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

    // Fluxes
    if (do_reflux) {
         if (current) 
         {
           for (int i = 0; i < BL_SPACEDIM ; i++) 
             current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, 1);
         }
         if (fine) 
         {
           for (int i = 0; i < BL_SPACEDIM ; i++) {
	         fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.,FluxRegister::ADD);
         }
      }
    }

    grav_vector.clear();

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0)
    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << "OOPS -- EFFECTIVE CFL AT THIS LEVEL " << level
                      << " IS " << courno << '\n';

        amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

}
#endif
