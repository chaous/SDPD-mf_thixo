/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(sdpd/taitwater/isothermal/mf,PairSDPDTaitwaterIsothermalMf)

#else

#ifndef LMP_PAIR_SDPD_TAITWATER_ISOTHERMAL_MF_H
#define LMP_PAIR_SDPD_TAITWATER_ISOTHERMAL_MF_H

#include "pair.h"
#ifdef USE_ZEST
#include <random>
#include "zest.hpp"
#endif

namespace LAMMPS_NS {

class PairSDPDTaitwaterIsothermalMf: public Pair {
 public:
  PairSDPDTaitwaterIsothermalMf (class LAMMPS *);
  virtual ~PairSDPDTaitwaterIsothermalMf();
  virtual void compute (int, int);
  void settings (int, char **);
  void coeff (int, char **);
  virtual double init_one (int, int);
  virtual void init_style();

 protected:
  
  //Modified by ASE, adding new variables for parameters ra, rb, A and F_int models (F) from Tartakovsky. Include a variable "fluid_type" to select newtonian (1) or thixotropic (2)
  double temperature, pb, rhoeq, ra, rb, A, fluid_type;
  int F;
  
  //Modified by ASE, adding new variables wall_index, a, b, alpha, f0 for thixotropic model
  double **wall_index, **a_thixo, **b_thixo, **alpha_thixo, **f0_thixo;
  double *rho0, *soundspeed, *B;
  double **cut;
  
  //Modified by ASE, adding new variable s_boo for the shape function condition (1 for repulsive-atractive and 0 for just repulsive) 
  double **sb;  
  
  //Modified by NMC, adding new variable for boundary interaction flag
  double **viscosity; 
  double **bulk_viscosity;
  double **slip;
  double **bccoord;
  double **S;

  void allocate ();

  int first;
  
  unsigned int seed;
#ifdef USE_ZEST
  std::mt19937_64 generator;
  Ziggurat<zest::StandardNormal,std::mt19937_64> gaussian;
#else
  class RanMars *random;
#endif
};

}

#endif
#endif
