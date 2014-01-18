#include "../text1r.h"
#include <unistd.h>

main(){
  double ttc=0.8;
  double p=15;
  double nu0=826320;
  double r=0.3;
  int itype=0;
  int fd = 6;
  int msglev = -3;
  double omega=1, omega_v=0.2, kr=1;
  int i;

  text1r_init_(&ttc, &p, &nu0, &r, &itype);
  text1r_pars_.lo=5;
//  text1r_set_vortex_twisted_(&omega, &kr);
  text1r_set_vortex_cluster_(&omega, &omega_v);
//  text1r_set_vortex_uniform_(&omega, &omega_v);
  text1r_minimize_(&msglev);
//  text1r_print("result.dat");
}

