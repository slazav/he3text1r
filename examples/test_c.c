#include "../text1r.h"
#include <unistd.h>
#include <string.h>

main(){
  double ttc = 0.8;     // Temperature, T/Tc
  double p   = 15;      // Pressure, bar
  double nu0 = 826320;  // Larmor frequency, Hz
  double r   = 0.3;     // Conteiner raduis, cm
  int itype  = 0;       // Initial conditions: 0 - normal, 1 - with 90deg peak...
  int n      = 100;     // Number of points
  int msglev = -3;      // Message level: -3 silent ...
  double omega   = 1;   // Rotation velocity, rad/s
  double omega_v = 0.2; // Rotation velocity of vortex cluster, rad/s
  double kr      = 1;   // Parameter for a twisted vortex profile
  char *fname = "result.dat"; // file for output

  // Initialize texture calculation:
  text1r_init_(&ttc, &p, &nu0, &r, &n, &itype);

  // Set vortex and velocity profile if needed:
  text1r_pars_.lo=5; // lambda/omega is not set by default
  text1r_set_vortex_cluster_(&omega, &omega_v);
//  text1r_set_vortex_uniform_(&omega, &omega_v);
//  text1r_set_vortex_twisted_(&omega, &kr);

  // Do minimization:
  text1r_minimize_(&msglev);

  // Print results to file:
  text1r_print_(fname, strlen(fname));
}
