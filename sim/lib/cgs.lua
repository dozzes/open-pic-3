-- cgs.lua -- Physical constants in CGS-Gaussian units
--
-- WHY CGS-Gaussian?
-- Plasma physics traditionally uses CGS-Gaussian units because Maxwell's
-- equations take a symmetric form and the speed of light c appears explicitly,
-- making electromagnetic scaling transparent.
-- In CGS: [E] = statV/cm, [B] = Gauss, [charge] = statcoulomb (esu).

return {
    c  = 2.9979e+10,   -- speed of light           [cm/s]
    e  = 4.8032e-10,   -- elementary charge (esu)  [statcoulomb]
    me = 9.1e-28,      -- electron mass             [g]
    mp = 1.6726e-24,   -- proton mass               [g]
}
