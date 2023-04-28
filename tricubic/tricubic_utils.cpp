#include "ltricubic_utils.h"

int ijk2n(int i, int j, int k) {
  return(i+4*j+16*k);
}

void point2xyz(int p, int *x, int *y, int *z) {
  switch (p) {
  case 0: *x=0; *y=0; *z=0; break;
  case 1: *x=1; *y=0; *z=0; break;
  case 2: *x=0; *y=1; *z=0; break;
  case 3: *x=1; *y=1; *z=0; break;
  case 4: *x=0; *y=0; *z=1; break;
  case 5: *x=1; *y=0; *z=1; break;
  case 6: *x=0; *y=1; *z=1; break;
  case 7: *x=1; *y=1; *z=1; break;
  default:*x=0; *y=0; *z=0;
  }
}
