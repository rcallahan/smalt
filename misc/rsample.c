/** Draw random samples from probability distributions */

#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include "elib.h"
#include "rsample.h"

int rsampleNormal(double *val, int n, double mu, double va)
{
  int i, n_steps = (n+1)/2;
  double ua, ub, lua, r, t;
  double PI = acos(-1.0);
  double v;

  /* uses the Box-Muller method */

  if (va < 0)
    return ERRCODE_ARGRANGE;
  v = sqrt(va);

  for (i=0; i<n_steps; i++) {
    ua = RANDRAW_UNIFORM_1();
    if (ua < DBL_EPSILON) ua = (double) 1;
    ub = RANDRAW_UNIFORM_1();
    if (ub < DBL_EPSILON) ub = (double) 1;
    lua = -2*log(ua);
    if (lua < 0) lua = 0;
    r = sqrt(lua);
    t = 2*PI*ub;
    *val++ = v*r*cos(t)+mu;
    n--;
    if (n<1) break;
    *val++ = v*r*sin(t)+mu;
    n--;
    if (n<1) break;
  }

  return (n>0)? ERRCODE_FAILURE: ERRCODE_SUCCESS;
}

int rsampleGeometric(int *val, int n, double p)
{
  int i;
  double rv, dr, dp, ldp;

  if (p < ((double) 0) || p >= ((double) 1))
    return ERRCODE_ARGRANGE;

  dp = 1.0 - p;
  if (dp < DBL_EPSILON) dp = DBL_EPSILON;
  ldp = log(dp);
  if (ldp > -DBL_EPSILON)
    ldp = - DBL_EPSILON;

  for (i=0; i<n; i++) {
    dr = RANDRAW_UNIFORM_1();
    if (dr < DBL_EPSILON) dr = (double) 1;
    rv = log(dr)/ldp;
    if (rv > INT_MAX)
      val[i] = INT_MAX;
    else
      val[i] = (int) rv;
  }
 
  return ERRCODE_SUCCESS;
}
