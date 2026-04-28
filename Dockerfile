# syntax=docker/dockerfile:1.6
#
# SDPD-mf_thixo: LAMMPS 21Nov2023 with the SDPD multiphase + thixotropy patches
# from BCAM-CFD/SDPD-mf_thixo, plus the missing atom_vec_sph.cpp fix.
#
# Build:
#   docker build -t sdpd-mf .
#
# Run (single rank, output mounted to host):
#   docker run --rm -v "$(pwd)/results:/work" sdpd-mf run-example 01_Young-laplace
#
# See README.md for more usage.

ARG LAMMPS_URL=https://download.lammps.org/tars/lammps-21Nov2023.tar.gz

# ---------- Stage 1: build LAMMPS with patches ----------
FROM ubuntu:24.04 AS builder

ARG LAMMPS_URL
ARG NPROC=4

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential \
        openmpi-bin libopenmpi-dev \
        wget ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

RUN wget -q "${LAMMPS_URL}" -O lammps.tar.gz && \
    tar -xzf lammps.tar.gz && \
    rm lammps.tar.gz && \
    mv lammps-* lammps

# Patched base files (from the upstream repo's lammps-21Nov23.tar.gz)
COPY atom.cpp atom.h atom_vec.cpp atom_vec.h /build/lammps/src/

# Patched SPH files: fix_sph + the atom_vec_sph fix that's missing upstream
COPY fix_sph.cpp fix_sph.h atom_vec_sph.cpp atom_vec_sph.h /build/lammps/src/SPH/

# New pair style
COPY pair_sdpd_taitwater_isothermal_mf.cpp pair_sdpd_taitwater_isothermal_mf.h \
     /build/lammps/src/DPD-SMOOTH/

# Register the new pair style with the DPD-SMOOTH package
RUN sed -i '/^action pair_sdpd_taitwater_isothermal\.cpp$/a\
action pair_sdpd_taitwater_isothermal_mf.h\
action pair_sdpd_taitwater_isothermal_mf.cpp' \
    /build/lammps/src/DPD-SMOOTH/Install.sh

# Build with MPI support
RUN cd /build/lammps/src && \
    make yes-SPH && \
    make yes-DPD-SMOOTH && \
    make -j"${NPROC}" mpi

# ---------- Stage 2: slim runtime image ----------
FROM ubuntu:24.04

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi3t64 && \
    rm -rf /var/lib/apt/lists/* && \
    useradd -m -u 1000 sdpd

COPY --from=builder /build/lammps/src/lmp_mpi /usr/local/bin/lmp_mpi
COPY --chown=sdpd:sdpd ["Numerical examples", "/opt/sdpd/examples"]
COPY --chown=sdpd:sdpd docker/run-example.sh /usr/local/bin/run-example
RUN chmod +x /usr/local/bin/run-example

USER sdpd
WORKDIR /work

# Default command lists available examples; override on the docker-run line.
CMD ["bash", "-c", "echo 'Available examples:'; ls /opt/sdpd/examples; echo; echo 'Run one with:'; echo '  docker run --rm -v \"$(pwd)/results:/work\" sdpd-mf run-example 01_Young-laplace'"]
