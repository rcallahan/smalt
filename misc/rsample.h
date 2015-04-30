/** Draw random samples from probability distributions */

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef RSAMPLE_H
#define RSAMPLE_H

#include "randef.h"

  int rsampleNormal(double *val, int n, double mu, double va);
  /**< Generate a sample of n N(mu,va) normal distributed values.
   * \param val Array of n elements returning the sample.
   * \param n Sample size (number of elements in array val).
   * \param mu Distribution mean.
   * \param va Distribution variance.
   */

  int rsampleGeometric(int *val, int n, double p);
  /**< Generate a sample of size n from the geometric distribution
   * \param val Array of n elements returning the sample.
   * \param  n Sample size (number of elements in array val).
   * \param p Probability 0 < p < 1.
   */

#endif
#ifdef __cplusplus
}
#endif
