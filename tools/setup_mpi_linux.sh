#!/usr/bin/env bash
# setup_mpi_linux.sh
# Installs an MPI implementation on Linux if not already present.
# Detects the package manager and installs OpenMPI or MPICH accordingly.
#
# Usage:
#   bash tools/setup_mpi_linux.sh
#
# Requires: sudo privileges.

set -euo pipefail

have() { command -v "$1" &>/dev/null; }

# -----------------------------------------------------------------------
# Check if MPI is already usable
# -----------------------------------------------------------------------

check_installed() {
    if have mpiexec && have mpicc; then
        echo "MPI already installed: $(mpiexec --version 2>&1 | head -1)"
        echo "  mpiexec : $(command -v mpiexec)"
        echo "  mpi.h   : $(find /usr/include /usr/local/include -name mpi.h 2>/dev/null | head -1)"
        return 0
    fi
    return 1
}

if check_installed; then
    echo ""
    echo "Nothing to do. Build with:"
    echo "  cmake -B build -DOPENPIC_ENABLE_MPI=ON"
    exit 0
fi

# -----------------------------------------------------------------------
# Detect package manager
# -----------------------------------------------------------------------

echo ""
echo "=== MPI setup for Linux ==="
echo ""

if have apt-get; then
    PKG_MGR="apt"
elif have dnf; then
    PKG_MGR="dnf"
elif have yum; then
    PKG_MGR="yum"
elif have zypper; then
    PKG_MGR="zypper"
elif have pacman; then
    PKG_MGR="pacman"
else
    echo "ERROR: No supported package manager found (apt/dnf/yum/zypper/pacman)."
    echo "Install OpenMPI or MPICH manually, then re-run."
    exit 1
fi

echo "Package manager: $PKG_MGR"
echo ""

# -----------------------------------------------------------------------
# Install
# -----------------------------------------------------------------------

case "$PKG_MGR" in
apt)
    sudo apt-get update -qq
    sudo apt-get install -y libopenmpi-dev openmpi-bin
    ;;
dnf|yum)
    sudo "$PKG_MGR" install -y openmpi-devel
    # OpenMPI uses environment modules on RHEL/Fedora
    MPI_MODULE=$(module avail mpi 2>&1 | grep -o 'mpi/openmpi[^ ]*' | head -1 || true)
    if [ -n "$MPI_MODULE" ]; then
        echo ""
        echo "NOTE: load the MPI module before building:"
        echo "  module load $MPI_MODULE"
    fi
    ;;
zypper)
    sudo zypper install -y openmpi4-devel
    ;;
pacman)
    sudo pacman -Sy --noconfirm openmpi
    ;;
esac

echo ""
echo "Verifying ..."
if check_installed; then
    echo ""
    echo "Done. Build with:"
    echo "  cmake -B build -DOPENPIC_ENABLE_MPI=ON"
    echo ""
    echo "Run with MPI:"
    echo "  mpiexec -n 4 ./build/bin/open-pic -mpi main.lua"
else
    echo ""
    echo "WARNING: mpiexec/mpicc not found in PATH after install."
    echo "You may need to reload your shell or load an environment module."
    echo "Then build with: cmake -B build -DOPENPIC_ENABLE_MPI=ON"
fi
