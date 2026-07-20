#pragma once

#include "grid/vector_3d.h"

#include <complex>
#include <vector>

namespace PIC {

// Scratch buffers for the field/electron solvers. Owned by simulate() and
// passed down by reference. Replaces function-local static vectors, which
// were never freed, blocked re-entrancy (two runs in one process shared the
// same buffers), and -- for the hyper-resistivity filter -- abused a
// DensityGrid's UP field as Laplacian scratch storage.
// The vectors keep their capacity between steps, so the per-step allocation
// behavior is identical to the old statics.
struct SolverWorkspace
{
    // PSTD spectral buffers (calc_magnetic_field_pstd_half_time):
    // real-space input, spectrum, dB spectrum, dB real-space.
    std::vector<std::complex<double>> pstd_Ex_r, pstd_Ey_r, pstd_Ez_r;
    std::vector<std::complex<double>> pstd_Ex_c, pstd_Ey_c, pstd_Ez_c;
    std::vector<std::complex<double>> pstd_dBx_c, pstd_dBy_c, pstd_dBz_c;
    std::vector<std::complex<double>> pstd_dBx_r, pstd_dBy_r, pstd_dBz_r;

    // Hyper-resistivity Laplacian of B (apply_hyper_resistivity),
    // flat [nx*ny*nz] layout, index (x*ny + y)*nz + z.
    std::vector<DblVector> lap_B;

    // Binomially smoothed Te fed into eta(Te) (update_Spitzer_coefficients)
    // and the ping-pong pass buffers of smooth_te_binomial.
    std::vector<double> te_smooth;
    std::vector<double> te_pass_a;
    std::vector<double> te_pass_b;
};

} // namespace PIC
