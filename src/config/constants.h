#pragma once

namespace PIC {

/*********************
 * Physical constants *
 *********************/
class Constants
{
  public:
    // Light velocity in vacuum (cm/sec)
    static constexpr double c() { return 2.99792458e+10; }

    // Electron charge in statcoulombs (Fr)
    static constexpr double e() { return 4.80320427e-10; }

    // Proton mass (g)
    static constexpr double mp() { return 1.6726219e-24; }

    // Real electron mass (g)
    // me ~ 9.109e-28, used for inertia in Ohm's law
    static constexpr double me() { return 9.10938356e-28; }

    static constexpr double pi() { return 3.14159265358979323846; }

    // --- Derived constants used in OpenPIC formulas ---

    // Proton charge/mass ratio: e/mp
    static constexpr double e_mp() { return 2.87175e+14; }

    // Coefficient for current density in Maxwell equations: c / (4 * pi * e)
    // Used to convert current/fields
    static constexpr double c_4pi_e() { return 4.96679926e+18; }

    // c / (4 * pi) — prefactor in Ampere's law: J = (c/4π) * curl(B)
    // Units: cm/s
    static constexpr double c_4pi() { return 2.38567258e+09; }

    // Electron charge to mass ratio (needed for electron inertia)
    static constexpr double e_me() { return 5.2728e+17; }

    // Boltzmann constant (erg/K)
    static constexpr double k_B() { return 1.380649e-16; }
};

} // namespace PIC
