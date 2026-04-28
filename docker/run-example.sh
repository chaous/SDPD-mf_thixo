#!/usr/bin/env bash
# Copy a numerical example into the (host-mounted) working directory and run it.
# Output dumps and log.lammps end up alongside the input on the host.
#
# Usage:
#   run-example <example_name> [extra lmp_mpi args...]
#
# Example:
#   run-example 01_Young-laplace
#   run-example 05_Poiseuille_one_phase

set -euo pipefail

EXAMPLES_DIR=/opt/sdpd/examples

if [ "$#" -lt 1 ]; then
    echo "Usage: run-example <example_name> [extra lmp_mpi args...]"
    echo
    echo "Available examples:"
    ls "${EXAMPLES_DIR}"
    exit 1
fi

NAME="$1"
shift

SRC="${EXAMPLES_DIR}/${NAME}"
if [ ! -d "${SRC}" ]; then
    echo "Example not found: ${NAME}"
    echo
    echo "Available examples:"
    ls "${EXAMPLES_DIR}"
    exit 1
fi

DST="/work/${NAME}"
mkdir -p "${DST}/dumps"
cp -r "${SRC}/." "${DST}/"
cd "${DST}"

INPUT=$(ls in.* 2>/dev/null | head -n 1)
if [ -z "${INPUT}" ]; then
    echo "No input file (in.*) found in ${NAME}" >&2
    exit 1
fi

echo "Running ${NAME}/${INPUT}"
echo "Output will be written to host volume mounted at /work/${NAME}"
echo
exec lmp_mpi -in "${INPUT}" "$@"
