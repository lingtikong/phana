#include <math.h>
#include "tricubic.h"
#include "ltricubic_utils.h"
#include "coeff.h"

char tricubic_version_stored[2048] = "0.1";

char *tricubic_version(void) {
  return(tricubic_version_stored);
}

void tricubic_pointID2xyz(int id, int *x, int *y, int *z) {
  point2xyz(id,x,y,z);
}

void tricubic_pointID2xyz(int id, double *x, double *y, double *z) {
  int x2,y2,z2;
  point2xyz(id,&x2,&y2,&z2);
  *x=(double)(x2);
  *y=(double)(y2);
  *z=(double)(z2);
}

void tricubic_get_coeff_stacked(double a[64], double x[64]) {
  int i,j;
  for (i=0;i<64;i++) {
    a[i]=(double)(0.0);
    for (j=0;j<64;j++) {
      a[i]+=A[i][j]*x[j];
    }
  }
}

void tricubic_get_coeff(double a[64], double f[8], double dfdx[8], double dfdy[8], double dfdz[8], double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8]) {
  int i;
  double x[64];
  for (i=0;i<8;i++) {
    x[0+i]=f[i];
    x[8+i]=dfdx[i];
    x[16+i]=dfdy[i];
    x[24+i]=dfdz[i];
    x[32+i]=d2fdxdy[i];
    x[40+i]=d2fdxdz[i];
    x[48+i]=d2fdydz[i];
    x[56+i]=d3fdxdydz[i];
  }
  tricubic_get_coeff_stacked(a,x);
}

double tricubic_eval(double a[64], double x, double y, double z) {
  int i,j,k;
  double ret=(double)(0.0);
  /* TRICUBIC EVAL
     This is the short version of tricubic_eval. It is used to compute
     the value of the function at a given point (x,y,z). To compute
     partial derivatives of f, use the full version with the extra args.
  */
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      for (k=0;k<4;k++) {
	ret+=a[ijk2n(i,j,k)]*pow(x,i)*pow(y,j)*pow(z,k);
      }
    }
  }
  return(ret);
}

double tricubic_eval(double a[64], double x, double y, double z, int derx, int dery, int derz) {
  int i,j,k;
  double ret=(double)(0.0);
  double cont;
  int w;
  /* TRICUBIC_EVAL 
     The full version takes 3 extra integers args that allows to evaluate
     any partial derivative of f at the point
     derx=dery=derz=0 => f
     derx=2 dery=derz=0 => d2f/dx2
     derx=dery=derz=1 =? d3f/dxdydz
     NOTICE that (derx>3)||(dery>3)||(derz>3) => returns 0.0
     this computes   \frac{\partial ^{derx+dery+derz} d}{\partial x ^{derx} \partial y ^{dery} \partial z ^{derz}}
  */
  for (i=derx;i<4;i++) {
    for (j=dery;j<4;j++) {
      for (k=derz;k<4;k++) {
	cont=a[ijk2n(i,j,k)]*pow(x,i-derx)*pow(y,j-dery)*pow(z,k-derz);
	for (w=0;w<derx;w++) {
	  cont*=(i-w);
	}
	for (w=0;w<dery;w++) {
	  cont*=(j-w);
	}
	for (w=0;w<derz;w++) {
	  cont*=(k-w);
	}
	ret+=cont;
      }
    }
  }
  return(ret);
}
