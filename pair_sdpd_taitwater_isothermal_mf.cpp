/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://lammps.sandia.gov/, Sandia National Laboratories
 Steve Plimpton, sjplimp@sandia.gov

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
      Morteza Jalalvand (IASBS)  jalalvand.m AT gmail.com

    references: Espanol and Revenga, Phys Rev E 67, 026705 (2003)
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------
	@autor:  Nicolas Moreno (BCAM), nmoreno@bcamath.org 

  Modified version of SDPD including:
  
  a) Explicitly bulk viscosity as presented by Espanol and Revenga, including background pressure, and generalizing for 2 and 3 dimensions.
	
  b) Update non-slip in flat wall: This update included the option for non-slip boundary condition in a flat
  wall as described in Bian and Ellero. The correction can be turned on and off using the variable slip[k][l]
  where the index correspond to type of particles k and l. Thus a fluid particle type k, can interact with non-slip bc
  with a wall particle type l.

  c) Update MF: This updated multiphase (mf) version is based on the description presented in
  Lei et al. 2016 (10.1103/PhysRevE.94.023304) including a pair-force contribution between particles.

  d) Transient viscosity model: This updated inclues a thixitropic (viscosity transient) model to simulate complex multiphase flows
-------------------------------------------------------------------------*/

#include "pair_sdpd_taitwater_isothermal_mf.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <math.h> 
#include <iostream>

#ifndef USE_ZEST
#include "random_mars.h"
#endif

using namespace LAMMPS_NS;

static const double sqrt_2_inv = std::sqrt(0.5);

/* ---------------------------------------------------------------------- */

PairSDPDTaitwaterIsothermalMf::PairSDPDTaitwaterIsothermalMf (LAMMPS *lmp)
: Pair (lmp) {
  restartinfo = 0;
  first = 1;
  single_enable =0;
}

/* ---------------------------------------------------------------------- */

PairSDPDTaitwaterIsothermalMf::~PairSDPDTaitwaterIsothermalMf () {
  if (allocated) {
    memory->destroy (setflag);
    memory->destroy (cutsq);

    memory->destroy (cut);
    memory->destroy (rho0);
    memory->destroy (soundspeed);
    memory->destroy (B);
    memory->destroy (viscosity);
    memory->destroy (bulk_viscosity);    
    memory->destroy (slip);
    memory->destroy (bccoord);

  }
}

/* ---------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalMf::compute (int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz,  evdwl,fpair;

  /* Interface force and shape factor to be used in the computation of multiphase. Following Lei2016 */
  double finter, phi, lambda, rasq, rbsq; 
  
  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc;
  double h, ih, ihsq, velx, vely, velz;
  double rsq, tmp, wfd, delVdotDelR;
  double prefactor, wiener[3][3], f_random[3];

  /* DNS added tmp variable for the force */
  double fxtmp, fytmp, fztmp;
  /* DNS end edit for ev_tally -> ev_tally_xyz */

  /* NMC added variable for thermal noises*/
  double a,b,Aij,Bij,trace,prefKernel,factor_sdpd;
  /* NMC end edit for ev_tally -> ev_tally_xyz */

  //NMC added distance from particle to flat boundary at fix y position
  // This is a test further implementation to any boundary needed
  double dbi, dbj, interslip;
  //double yboundary = bccoord;// domain->boxlo[1];

  evdwl = 0.0;
  ev_init(eflag, vflag);
  
  /* ASE added termporal variables for thixotropic model*/
  double gamma,visco,bulk_visco,eta_12,bulk_eta_12;
  double lambda_time_scale,beta_thixo,exp_thixo,f_a,f_b,f_t,y_s;

  /* ASE added termporal variable for viscosity*/
  double visco_temp, bulk_visco_temp, gammadot_tmp;
  
  double *f_dst;   // parámetro de orden φᵢ
  double df_dst;

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  //add by NMC for enabling lj special
  double *special_lj = force->special_lj;

  int dimension = domain->dimension;

  double **gradv = atom->gradv;
  double *gammadot = atom->gammadot;
  double *fintx = atom->fintx;
  double *finty = atom->finty;  

  double curTime = update->dt * update->ntimestep;
  double dtinv = 1.0 / update->dt;
  double kBoltzmann = force->boltz;
  double dim;// = (double)dimension;	

if(dimension==3){
  dim=3.0;
} else {
  dim=2.0;
}


  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


//NMC distances used by the pair-wise potential for MF 
rbsq = rb*rb;
rasq = ra*ra;

// prefactor for kernel in terms of the box dimensions
prefKernel = -25.066903536973515383e0*(dim-2.0) - 19.098593171027440292e0*(3.0-dim); //first term for 3d second 2d

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    // compute pressure of atom i with Tait EOS
    tmp = rho[i] / rho0[itype];
    fi = tmp * tmp * tmp;
    //Modified by NMC including background pressure
    fi = (B[itype] * (fi * fi * tmp - 1.0)+pb) / (rho[i] * rho[i]);

    
    //Modified by ASE including gamma_dot declaration	
    gammadot_tmp = gammadot[i];  


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_sdpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      dbj = (bccoord[itype][jtype] - x[j][1]);
      dbi = (bccoord[itype][jtype] - ytmp);


      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;

        double r = sqrt (rsq);

        // ASE adds the first three F_int models (F) defined by Tartakovsky_2016
        
        if (F==1){
        phi  = cos((3*3.14159265359)/(2*h)*r)/r;     			// Shape function for pairwise force F_int_1 
        finter = -(-S[itype][jtype]*phi);		        // Pairwise force F_int_1          
        } else if (F==2){
        phi  = (-A*exp(-rsq/(2.0*rasq))+exp(-rsq/(2.0*rbsq))*(1.0-sb[itype][jtype]))/r;      // Shape function for pairwise force F_int_2 
        finter = -S[itype][jtype]*phi;		        // Pairwise force F_int_2          
        } else if (F==3){
        phi  = -A*exp(-rsq/(2.0*rasq)) + exp(-rsq/(2.0*rbsq))*(1.0-sb[itype][jtype]);    // Shape function for pairwise force F_int_3
        finter = -S[itype][jtype]*phi;		        // Pairwise force F_int_3                   
        }
        
        
        finter*=factor_sdpd;
	
	// Interaction forces due to the new pairwise               
        fintx[i] = finter*delx;
        finty[i] = finter*dely;

	//Kernel definition
        wfd = h - r;
        wfd = prefKernel * wfd * wfd * ihsq * ihsq * ihsq * (1*(3.0-dim)+ih*(dim-2.0));
        
        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
        fj = (B[jtype] * (fj * fj * tmp - 1.0)+pb) / (rho[j] * rho[j]);


       // Non-slip in flat wall
          interslip = fabs(dbj/dbi);
          if (interslip>0.5 || isnan(interslip)){ interslip = 0.5;}

         velx = vxtmp - (1.0-slip[itype][jtype])*v[j][0] + slip[itype][jtype]*((interslip)*(vxtmp-v[j][0])-v[j][0]);
         vely = vytmp - (1.0-slip[itype][jtype])*v[j][1] + slip[itype][jtype]*((interslip)*(vytmp-v[j][1])-v[j][1]);
         velz = vztmp - (1.0-slip[itype][jtype])*v[j][2] + slip[itype][jtype]*((interslip)*(vztmp-v[j][2])-v[j][2]);


        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;


	//ASE add the calculation of the U gradient

        gradv[i][0] += -(jmass  * wfd / (rho[j])) * velx*delx;
        gradv[i][1] += -(jmass  * wfd / (rho[j])) * vely*dely;
        gradv[i][2] += -(jmass  * wfd / (rho[j])) * velz*delz;
        gradv[i][3] += -(jmass  * wfd / (rho[j])) * velx*dely;
        gradv[i][4] += -(jmass  * wfd / (rho[j])) * velx*delz;
        gradv[i][5] += -(jmass  * wfd / (rho[j])) * vely*delz;
        gradv[i][6] += -(jmass  * wfd / (rho[j])) * vely*delx;
        gradv[i][7] += -(jmass  * wfd / (rho[j])) * velz*delx;
        gradv[i][8] += -(jmass  * wfd / (rho[j])) * velz*dely;


//------Thixotropic model--------------------------------------------------------------------------------------

//Boolen operator to compute viscosity (0 if i=j and 1 if i =/ j)
        double f12 = (i == j) ? 0 : 1;


//  ASE added differents models for viscosity

// fluid_type==1 for Newtonian fluid
        if (fluid_type==1){
        
        	eta_12= 2*(viscosity[itype][jtype]*viscosity[itype][jtype]) / (viscosity[itype][jtype] + viscosity[itype][jtype]);
		bulk_eta_12 = 0.999*eta_12*(2.0*dim-1.0)/dim;
	
		visco = viscosity[itype][jtype]*(1.0-f12) + eta_12*(f12);
		bulk_visco = (0.999*visco*(2.0*dim-1.0)/dim)*(1.0-f12) + bulk_eta_12*(f12);

// fluid_type==2 for Thixotropic model
        } else if (fluid_type==2){


  	beta_thixo = (b_thixo[itype][jtype]/a_thixo[itype][jtype]);
                lambda_time_scale = (1.0/a_thixo[itype][jtype]);
                exp_thixo = exp(-(1+beta_thixo*gammadot_tmp)*curTime/lambda_time_scale);

                f_a = (1.0/(1.0+beta_thixo*gammadot_tmp))*(1.0 - exp_thixo);
                f_b = f0_thixo[itype][jtype]*exp_thixo;
                f_t = f_a + f_b;

//Viscosity calculation for i =/ j
                eta_12= 2*((viscosity[itype][jtype]*(1.0+alpha_thixo[itype][jtype]*f_t))*(viscosity[itype][jtype]*(1.0+alpha_thixo[itype][jtype]*f_t))) / ((viscosity[itype][jtype]*(1.0+alpha_thixo[itype][jtype]*f_t))+(viscosity[itype][jtype]*(1.0+alpha_thixo[itype][jtype]*f_t)));

                bulk_eta_12 = 0.999*eta_12*(2.0*dim-1.0)/dim;

                visco = (viscosity[itype][jtype]*(1.0+alpha_thixo[itype][jtype]*f_t))*(1.0-f12) + eta_12*(f12);
                bulk_visco = (0.999*visco*(2.0*dim-1.0)/dim)*(1.0-f12) + bulk_eta_12*(f12);
        } 	

//-----------------------------------------------------------------------------------------------
	
//NMC added amplitude of thermal noise, accounting for different dimensionality         // Espanol Viscosity (Espanol, 2003)
        a = (2.0-1.0/dim)*(visco) - bulk_visco; 
        b = (2.0+dim)/dim * (visco)+(2.0+dim)*bulk_visco - a*(2.0*dim-4.0)/(2.0*dim); //b+1/3a  eqs 61 Espanol2003, //for 2D a contribution on b vanishes

        fvisc = imass * jmass * wfd / (rho[i]*rho[j]);
        fvisc*=factor_sdpd;


        // total pair force
        fpair = -imass * jmass * (fi + fj) * wfd;
        fpair*=factor_sdpd;


 // random force calculation
        // independent increments of a Wiener process matrix


        // NMC redefine matrix as eqs 50 and 51 
/*
#ifdef USE_ZEST
        wiener[0][0] = gaussian (generator);
        wiener[1][1] = gaussian (generator);
        wiener[2][2] = gaussian (generator);

        wiener[0][1] = wiener[1][0] = 0.5 * (gaussian (generator)+gaussian (generator));
        wiener[0][2] = wiener[2][0] = 0.5 * (gaussian (generator)+gaussian (generator));
        wiener[1][2] = wiener[2][1] = 0.5 * (gaussian (generator)+gaussian (generator));
#else
*/
        wiener[0][0] = random->gaussian ();
        wiener[1][1] = random->gaussian ();

        wiener[0][1] = wiener[1][0] = random->gaussian ();//0.5 * (random->gaussian ()+random->gaussian ());
        wiener[0][2] = wiener[2][0] = random->gaussian ();//0.5 * (random->gaussian ()+random->gaussian ());
      //  wiener[1][2] = wiener[2][1] = random->gaussian ();//0.5 * (random->gaussian ()+random->gaussian ());

//#endif

        /* NMC random tensor for general case with kin and bulk viscosities*/
        prefactor = -4.0 * kBoltzmann*temperature * fvisc * dtinv;
        Aij = sqrt(prefactor * a)/r;
        Bij = sqrt(prefactor * b * dim/2.0)/r;

        if (dimension==3)
        {
        	 wiener[2][2] = random->gaussian ();
        	 wiener[1][2] = wiener[2][1] = random->gaussian ();//0.5 * (random->gaussian ()+random->gaussian ());
             wiener[2][2] = random->gaussian ();

          trace = (wiener[0][0]+wiener[1][1]+wiener[2][2])/3.0;

        wiener[0][0] = wiener[0][0]-trace;
        wiener[1][1] = wiener[1][1]-trace;
        wiener[2][2] = wiener[2][2]-trace;

        /* NMC random tensor if bulk viscosity*/
        f_random[0] = factor_sdpd*(Aij*(wiener[0][0]*delx + wiener[0][1]*dely + wiener[0][2]*delz) + Bij*trace*delx);
        f_random[1] = factor_sdpd*(Aij*(wiener[1][0]*delx + wiener[1][1]*dely + wiener[1][2]*delz) + Bij*trace*dely);
        f_random[2] = factor_sdpd*(Aij*(wiener[2][0]*delx + wiener[2][1]*dely + wiener[2][2]*delz) + Bij*trace*delz);

        // NMC including separate a and b terms eq 61 Espanol                                                           ///Multiphase
        f[i][0] += delx * fpair + (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc + f_random[0] + finter*delx;
        f[i][1] += dely * fpair + (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc + f_random[1] + finter*dely;
        f[i][2] += delz * fpair + (a*velz + (b+a*(2.0*dim-4.0)/(2.0*dim))*delz * delVdotDelR / rsq) * fvisc + f_random[2] + finter*delz;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        if (newton_pair || j < nlocal) {
          f[j][0] -= (delx * fpair + (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc + f_random[0]) + finter*delx;
          f[j][1] -= (dely * fpair + (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc + f_random[1]) + finter*dely;
          f[j][2] -= (delz * fpair + (a*velz + (b+a*(2.0*dim-4.0)/(2.0*dim))*delz * delVdotDelR / rsq) * fvisc + f_random[2]) + finter*delz;
          drho[j] += imass * delVdotDelR * wfd;
        }
          /* DNS edit for temp forces for ev_tally_xyz*/
        fxtmp = delx*fpair + (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc + f_random[0] + finter*delx;
        fytmp = dely*fpair + (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc + f_random[1] + finter*dely;
        fztmp = delz*fpair + (a*velz + (b+a*(2.0*dim-4.0)/(2.0*dim))*delz * delVdotDelR / rsq) * fvisc + f_random[2] + finter*delz; 
        } 
        else {

         trace = (wiener[0][0]+wiener[1][1])/2.0;
         wiener[0][0] = wiener[0][0]-trace;
         wiener[1][1] = wiener[1][1]-trace;

        /* NMC random tensor if bulk viscosity*/
        f_random[0] = factor_sdpd*(Aij*(wiener[0][0]*delx + wiener[0][1]*dely) + Bij*trace*delx);
        f_random[1] = factor_sdpd*(Aij*(wiener[1][0]*delx + wiener[1][1]*dely) + Bij*trace*dely);
        f_random[2] = 0;

        f[i][0] += delx * fpair + (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc + f_random[0] + finter*delx;
        f[i][1] += dely * fpair + (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc + f_random[1] + finter*dely;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;
        /*NMC edit forces sign was not properly pass the parenthesis was added  */
        if (newton_pair || j < nlocal) {
          f[j][0] -= (delx * fpair + (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc + f_random[0]) + finter*delx;
          f[j][1] -= (dely * fpair + (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc + f_random[1]) + finter*dely;
          
          drho[j] += imass * delVdotDelR * wfd;
          
          
          //Approximating gradient of the velocity at that point using the extrapolated velocity vest
          //needed to proper stimation of gradient with ghost atoms
          gradv[j][0] += -(imass * wfd / (rho[i])) * delx*velx;
          gradv[j][1] += -(imass * wfd / (rho[i])) * dely*vely;
          gradv[j][2] += -(imass * wfd / (rho[i])) * delz*velz;
          gradv[j][3] += -(imass * wfd / (rho[i])) * velx*dely;
          gradv[j][4] += -(imass * wfd / (rho[i])) * velx*delz;
          gradv[j][5] += -(imass * wfd / (rho[i])) * vely*delz;
          gradv[j][6] += -(imass * wfd / (rho[i])) * vely*delx; 
          gradv[j][7] += -(imass * wfd / (rho[i])) * velz*delx;
          gradv[j][8] += -(imass * wfd / (rho[i])) * velz*dely;
        }

        /* NMC edit for temp forces for ev_tally_xyz. delR factor out*/
        fxtmp = delx*fpair + (a*velx + (b+a*(2.0*dim-4.0)/(2.0*dim))*delx * delVdotDelR / rsq) * fvisc + f_random[0] + finter*delx;
        fytmp = dely*fpair + (a*vely + (b+a*(2.0*dim-4.0)/(2.0*dim))*dely * delVdotDelR / rsq) * fvisc + f_random[1] + finter*dely;
        fztmp = 0.0;
        } 


        /* NMC including the computation of energy term due pairs - only conservative*/
        if (eflag)
          evdwl = 0.5*(imass * jmass * (fi + fj) * cut[itype][jtype] * wfd*wfd);
          evdwl *= 1.;

        if (evflag)
          ev_tally_xyz (i, j, nlocal, newton_pair, evdwl, 0.0, fxtmp, fytmp, fztmp, delx, dely, delz);
          //ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute ();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalMf::allocate () {
  allocated = 1;
  int n = atom->ntypes;

  memory->create (setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create (cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create (rho0, n + 1, "pair:rho0");
  memory->create (soundspeed, n + 1, "pair:soundspeed");
  memory->create (B, n + 1, "pair:B");  
  memory->create (viscosity, n+1, n + 1, "pair:viscosity");
  memory->create (bulk_viscosity, n+1, n + 1, "pair:bulk_viscosity");
  memory->create (a_thixo, n+1, n + 1, "pair:a_thixo");
  memory->create (wall_index, n+1, n + 1, "pair:wall_index");  
  memory->create (b_thixo, n+1, n + 1, "pair:b_thixo");
  memory->create (alpha_thixo, n+1, n + 1, "pair:alpha_thixo");
  memory->create (f0_thixo, n+1, n + 1, "pair:f0_thixo");       
  memory->create (cut, n + 1, n + 1, "pair:cut");
  memory->create (slip, n + 1, n + 1, "pair:slip");
  memory->create (bccoord, n+1, n + 1, "pair:bccoord");
  memory->create (S, n+1, n + 1, "pair:S");
  memory->create (sb, n+1, n + 1, "pair:sb");

}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalMf::settings (int narg, char **arg) {
  if (narg != 8 && narg != 9)
    error->all (FLERR, "Illegal number of arguments for "
                "pair_style sdpd/taitwater/isothermal");

  temperature = utils::numeric(FLERR, arg[0], false, lmp);
  //viscosity = utils::numeric(FLERR, arg[1], false, lmp);
  /* NMC adding bulk viscosity and background pressure as variable*/
  //bulkViscosity = utils::numeric (FLERR, arg[2],false, lmp);
  pb = utils::numeric (FLERR, arg[1], false, lmp);
  rhoeq = utils::numeric(FLERR, arg[2],false,lmp);
  /* ASE adding factors ra,rb, A and F_int models (F) as variables */
  ra = utils::numeric(FLERR, arg[3],false,lmp);
  rb = utils::numeric(FLERR, arg[4],false,lmp);
  A = utils::numeric(FLERR, arg[5],false,lmp);
  F = utils::numeric(FLERR, arg[6],false,lmp);
  fluid_type = utils::numeric(FLERR, arg[7],false,lmp);

  if (temperature <= 0) error->all (FLERR, "Temperature must be positive");
  //if (viscosity < 0) error->all (FLERR, "Viscosity must be positive");
  //if (bulkViscosity < 0) error->all (FLERR, "Bulk viscosity must be positive");
  if (pb < 0) error->all (FLERR, "Background pressure must be positive");
  if (rhoeq < 0) error->all (FLERR, "Interfacial tension must be possitive");
  /* ASE adding conditional for model number F. It had to be an integer*/
  if (F < 1 || F > 3) error->all(FLERR, "Pairwise forces model must be an integer between 1 and 3");   
  /* ASE adding conditional for model number F. It had to be an integer*/
  if (fluid_type < 1 || fluid_type > 2) error->all(FLERR, "Viscosity model must be an integer between 1 and 2");   


  // seed is immune to underflow/overflow because it is unsigned
  seed = comm->nprocs + comm->me + atom->nlocal;
  if (narg == 9) seed += utils::inumeric(FLERR, arg[8], false, lmp);
#ifdef USE_ZEST
  generator.seed (seed);
#else
  random = new RanMars (lmp, seed);
#endif
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalMf::coeff (int narg, char **arg) {
  if (narg != 13 && narg != 16)
    error->all (FLERR, "Incorrect args for pair_style "
                "sph/taitwater/isothermal coefficients");

  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double rho0_one = utils::numeric(FLERR,arg[2], false, lmp);
  double soundspeed_one = utils::numeric(FLERR,arg[3], false, lmp);
  double cut_one = utils::numeric(FLERR,arg[4], false, lmp);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;
  double viscosity_one = utils::numeric(FLERR,arg[5], false, lmp);
  double bulk_viscosity_one = utils::numeric(FLERR,arg[6], false, lmp);
  /* ASE adding factors wall_index, a_thixo, b_thixo and alpha_thixo for the thixotropic model as variables*/
  double wall_index_one = utils::numeric(FLERR,arg[7], false, lmp);     
  double a_thixo_one = utils::numeric(FLERR,arg[8], false, lmp);   
  double b_thixo_one = utils::numeric(FLERR,arg[9], false, lmp);
  double alpha_thixo_one = utils::numeric(FLERR,arg[10], false, lmp);
  double f0_thixo_one = utils::numeric(FLERR,arg[11], false, lmp);
  
  double bccoord_one = 0.0; //So far only flat walls perpendicular to the y axis
  double slip_one = 0.0; //Default not slip interpolation between particles i j 
  double S_one = utils::numeric(FLERR,arg[12], false, lmp); //Default no interaction between same species 
  double sb_one = utils::numeric(FLERR,arg[13], false, lmp);

  if (narg==16){
  	///ADDITIONAL ERROR CATCHING LINES NEEDED TO ACCOUNT FOR PROPER POSITION OF THE WALL
    slip_one = utils::numeric(FLERR,arg[14], false, lmp);
    bccoord_one = utils::numeric(FLERR,arg[15],false,lmp);
   }
  if (rho0_one <= 0) error->all (FLERR, "Density must be positive");
  if (soundspeed_one <= 0) error->all (FLERR, "Sound speed must be positive");
  if (cut_one <= 0) error->all (FLERR, "Cutoff must be positive");
  if (a_thixo_one <= 0) error->all(FLERR, "Parameter a_thixo must be positive and different than 0");
  if (b_thixo_one < 0) error->all(FLERR, "Parameter b_thixo must be positive");
  //if (sb_one != 0.0 && sb_one != 1.0) error->all(FLERR, "Parameter sb must be either 0 or 1"); 
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      slip[i][j] = slip_one;
      bccoord[i][j] = bccoord_one;
      viscosity[i][j] = viscosity_one;
      bulk_viscosity[i][j] = bulk_viscosity_one;
      wall_index[i][j] = wall_index_one;
      a_thixo[i][j]=a_thixo_one;
      b_thixo[i][j]=b_thixo_one;
      alpha_thixo[i][j]=alpha_thixo_one;            
      f0_thixo[i][j]=f0_thixo_one;
      S[i][j] = S_one;
      sb[i][j] = sb_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init specific to this pair style
------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermalMf::init_style()
{
  if ((!atom->rho_flag) || (atom->drho == nullptr))
    error->all(FLERR,"Pair style dpd/taitwater/isothermal requires atom "
               "attributes rho and drho");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSDPDTaitwaterIsothermalMf::init_one (int i, int j) {
  if (setflag[i][j] == 0)
    error->all(FLERR,"Not all pair sdpd/taitwater/isothermal coeffs are set");

  cut[j][i] = cut[i][j];
  slip[j][i] = slip[i][j];
  bccoord[j][i] = bccoord[i][j];
  viscosity[j][i] = viscosity[i][j];
  bulk_viscosity[j][i] = bulk_viscosity[i][j];
  wall_index[j][i] = wall_index[i][j];      
  a_thixo[j][i] = a_thixo[i][j];  
  b_thixo[j][i] = b_thixo[i][j];
  alpha_thixo[j][i] = alpha_thixo[i][j];
  f0_thixo[j][i] = f0_thixo[i][j];

  S[j][i] = S[i][j];
  sb[j][i] = sb[i][j];

  return cut[i][j];
}

