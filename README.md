# SDPD-mf_thixo

![Image](image.png)

Smoothed Dissipative Particle Dynamics (SDPD) model for multiphase and thixotropic fluids implemented in LAMMPS.

# Description

An extended version of the standard SDPD model in LAMMPS (Jalalvand et al. 2021),  including:
 	
a) Update non-slip in flat wall: This update included the option for a non-slip boundary condition in a flat wall as described in Bian and Ellero 2012. The correction can be turned on and off using the variable slip[k][l], where the index corresponds to the type of particles k and l. Thus, a fluid particle type k can interact with a non-slip bc with a wall particle type l.

b) Update MF: This updated multiphase (mf) version is based on the description presented in Lei et al. 2016 (10.1103/PhysRevE.94.023304), including a pair-force contribution between particles.

c) Transient viscosity model: This update includes a thixotropic (viscosity transient) model to simulate complex multiphase flows

# Software version

This SDPD extension must be compiled in the 21-Nov-2023 release of the open source LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator).

# Quick start with Docker

If you prefer not to compile LAMMPS by hand, the included `Dockerfile` builds a self-contained image with the patched LAMMPS binary and all numerical examples. It works on any machine with Docker installed (Linux, macOS, Windows via Docker Desktop or WSL2).

## 1. Build the image

```
docker build -t sdpd-mf .
```

This downloads LAMMPS 21Nov2023 from the official mirror, applies all patches in this repository (including the `atom_vec_sph.cpp` fix described in [`atom_vec_sph.patch`](atom_vec_sph.patch)), and compiles `lmp_mpi` with MPI support. Build time is typically 5–15 minutes.

To use more cores during compilation:

```
docker build --build-arg NPROC=8 -t sdpd-mf .
```

## 2. Run a numerical example

Output (dump files, log) is written to a host directory you mount into the container.

```
mkdir -p results
docker run --rm -v "$(pwd)/results:/work" sdpd-mf run-example 01_Young-laplace
```

When the run finishes, you'll find `results/01_Young-laplace/dumps/*.lammpstrj` and `results/01_Young-laplace/log.lammps` on the host. Drag the `dumps/` folder into [OVITO](https://www.ovito.org/) or ParaView to visualise the trajectory.

## 3. List available examples

```
docker run --rm sdpd-mf
```

(The default command prints the list of cases bundled in `/opt/sdpd/examples`.)

## 4. Drop into a shell inside the container

```
docker run --rm -it -v "$(pwd)/results:/work" sdpd-mf bash
```

Useful for inspecting input files, running custom invocations of `lmp_mpi`, or debugging.

## 5. Notes

- The image runs as a non-root `sdpd` user (UID 1000); make sure the host directory you mount is writable by UID 1000, or run `docker run --user $(id -u):$(id -g) ...`.
- Multi-rank MPI runs (`mpirun -np N lmp_mpi …`) are not yet recommended — the `fields_comm` / `fields_border` lists for the new per-atom fields (`gammadot`, `gradv`, `fintx`, `finty`, `stmic`, `stfluid`) need confirmation from the authors before parallel correctness can be guaranteed. Stick to single-rank for now.

---

# Installation and Build flags (manual)

The installation process and build flags are detailed step by step below.

1. Download and unzip the .tar.gz file “lammps-21Nov23.tar.gz” found in this repository.
2. Replace the following files with the files from the .tar

   *atom.cpp*
   
   *atom.h*
   
   *atom_vec.cpp*
   
   *atom_vec.h*
   
   *fix_sph.cpp*
   
   *fix_sph.h*

4. Copy the files *pair_sdpd_taitwater_isothermal_mf.cpp* and *pair_sdpd_taitwater_isothermal_mf.h* in the folder "DPD-SMOOTH"
5. In the same folder, locate the “Install.sh” file. Include these two lines to indicate that you are going to install a new dependency:

   *action pair_sdpd_taitwater_isothermal_mf.h*
   
   *action pair_sdpd_taitwater_isothermal_mf.cpp*
   
6. LAMMPS must be compiled with the SPH and SDPD package enabled. To activate it, located in the “src” folder type:

   *make yes-SPH*
   
   *make yes-DPD-SMOOTH*
 
8. LAMMPS must be compiled using MPI support. Located in the “src” folder type:

   *make pu*
   
   *make mpi*

10. Run a numerical example to verify the installation. Copy any example to your local folder and, once located in the case folder, type:

   *lmp_mpi <in.sdpd_phase.2d*

# Numerical examples

In the folder called “Numerical examples” you will find a collection of all the numerical cases covered in this research. Detailed information on each case and how to run them can be found in the README file in each example folder. The complete list of cases studied is detailed below: 

1. Surface tension of a droplet - Used for static validation - See Figure 1_(a)
2. Retraction of a stretched droplet - Used for static validation - See Figure 1_(b)
3. Static contact angle between droplet and solid wall - Used for static validation - See Figure 1_(c)
4. Triple contact angle between three diferents droplets - Used for static validation - See detail in Figure 1_(c)
5. Poiseuille flow for one phase - Used for dynamic validation - See Figure A.14_(a) and A.14_(d)
6. Reverse Poiseuille flow for one phase - Used for dynamic validation - See Figure A.14_(b) and A.14_(e)
7. Flow around a cylinder - Used for dynamic validation - See Figure A.14_(c) and A.14_(f)
8. Poiseuille flow for two phases - Used for dynamic validation - See Figure 2_(a)
9. Taylor deformation vs capillary number - Used for dynamic validation - See Figure 2_(b)
10. Droplet break-up - Used for dynamic validation - See Figure 2_(c)
11. Thixotropic shear flow in a channel - Used for thixotropic validation - See Figure 3_(a) and 3_(b)
12. Liquid-Liquid Phase Separation (LLPS) - First exploratory case - See Figure 6 and 7
13. Poiseuille flow for two phases with one thixotropic phase - Second exploratory case - See Figure 8
14. Emulsion with thixotropic continuous phase and newtonian droplets - Second exploratory case - See Figure 9
15. Emulsion with newtonian continuous phase and thixotropic droplets - Second exploratory case - See Figure 10
16. Droplet dynamics in a periodically constricted channel - Third exploratory case - See Figure 11 and 12
17. Droplet merging using micro-devices - fourth exploratory case - See Figure 13
