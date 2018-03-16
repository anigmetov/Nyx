
#include "Nyx.H"
#include "Nyx_F.H"

using namespace amrex;
using std::string;

void
Nyx::sdc_reaction_step (Real time, Real dt, MultiFab& S_old, MultiFab& D_old, MultiFab& S_src, int sdc_iter)
{
    BL_PROFILE("Nyx::sdc_first_step()");
    Real half_dt = 0.5*dt;

    const Real a = get_comoving_a(time);
    const Real* dx = geom.CellSize();


        if(sdc_iter>1e-2)
          compute_new_temp();
    
#ifndef FORCING
    {
      const Real z = 1.0/a - 1.0;
      fort_interp_to_this_z(&z);
    }
#endif

    
#ifdef _OPENMP
#pragma omp parallel
#endif
    //    D_old.setVal(0,Sfnr_comp);
    //    D_old.setVal(-2,Ssnr_comp);
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        // Note that this "bx" is only the valid region (unlike for Strang)
      const Box& bx = mfi.tilebox();

        int  min_iter = 100000;
        int  max_iter =      0;

        integrate_state
                (bx.loVect(), bx.hiVect(), 
                 BL_TO_FORTRAN(S_old[mfi]),
                 BL_TO_FORTRAN(D_old[mfi]),
		 BL_TO_FORTRAN(S_src[mfi]),
                 dx, &time, &a, &dt, &min_iter, &max_iter, &sdc_split);

#ifndef NDEBUG
        if (S_old[mfi].contains_nan())
            amrex::Abort("state has NaNs after the first strang call");
#endif

    }
    
}
