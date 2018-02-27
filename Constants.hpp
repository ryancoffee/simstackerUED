#ifndef _CONSTANTS_H
#define _CONSTANTS_H

#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>

// setting up general constants ...as inline replacements via preprocessor //

//#define pi M_PI // 3.141592653589793
#define Eh 27.2113845 // eV/hartree
#define a0 0.5291772108 // AA/bohr
#define nmPau 0.05291772108 // nm/bohr
#define icm 8065.54445 // icm/eV
#define hc 1239.84190604789 // 197.326968*2*pi eV nm
#define C 2.99792458e10 // cm/s
#define C_nmPfs (GSL_CONST_MKSA_SPEED_OF_LIGHT * 1e-6 ) // C [nm/fs]
#define MC2 931.494e6 // eV
#define hbar 6.58211915e-16
#define kb 8.617343e-5 // eV / K
#define fsPau .02418884326505 // fs / au
#define icmPau (icm * Eh) // icm/hartree
#define amuPau (9.1093826e-31/1.66053886e-27) // unified AMU / au_mass
#define e2P4pieps0 14.3996445 //eV Ang
#define auPe2P4pieps0 (e2P4pieps0 /Eh/a0)  // hartree borh
#define aufor10PW 0.5336 // Atomic unit of field for 10^16 W cm^-1
#define auenergy 183.631526 // 4pi eps_0 * E_au^2 in units of  eV ang^-3
#define mbarn 1e-24 // Mb/cm^2
#define auIntensity 3.55e16 // W/cm^2
#define aupolarizability (a0*a0*a0) // ang^3

#endif

/*
Upinau(x,y) = 8*M_PI/137 * x/4/(y**2) // all in atomic units 
 = 0.1835 * x/ 4/y
 = 200 * x so for 3micron light, 1au of Up needs 0.005 au of intensity 
*/

