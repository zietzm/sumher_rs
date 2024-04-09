/*
Copyright 2022 Doug Speed.

    LDAK is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

    LDAK is distributed in the hope that they will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the fGNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with
LDAK.  If not, see <http://www.gnu.org/licenses/>.

*/

///////////////////////////

// summary functions

//////////////////////////
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "decomp.h"
#include "fortran_functions.h"
#include "sumfuns.h"

int solve_sums(double *stats, double *likes, double *cohers, double *influs,
               int num_parts, int gcon, int cept, int num_blocks, int length,
               int ncv, int *cvindex, double *cvexps, double *stags,
               double **svars, double **ssums, double *snss, double *schis,
               double tol, int maxiter, int chisol, int sflag)
// sflag=0 - normal, sflag=1 - first pass, sflag=2 - second pass
// sflag=3 - just get expectations and likelihood, sflag=4 - LDSC, sflag=5 -
// divide+updating
{
  // Print the first few of each input variable
  // printf("num_parts: %d\n", num_parts);
  // printf("gcon: %d\n", gcon);
  // printf("cept: %d\n", cept);
  // printf("num_blocks: %d\n", num_blocks);
  // printf("length: %d\n", length);
  // printf("ncv: %d\n", ncv);
  //
  // printf("cvindex: ");
  // for (int i = 0; i < ncv; i++) {
  //   printf("%d ", cvindex[i]);
  // }
  // printf("\n");
  //
  // printf("cvexps: ");
  // for (int i = 0; i < ncv; i++) {
  //   printf("%f ", cvexps[i]);
  // }
  // printf("\n");
  //
  // printf("stags: ");
  // for (int i = 0; i < 3; i++) {
  //   printf("%f ", stags[i]);
  // }
  // printf("\n");
  //
  // printf("svars: ");
  // for (int i = 0; i < num_parts; i++) {
  //   for (int k = 0; k < 3; k++) {
  //     printf("%f ", svars[i][k]);
  //   }
  // }
  // printf("\n");
  //
  // printf("ssums: ");
  // for (int i = 0; i < num_parts; i++) {
  //   for (int k = 0; k < (num_parts + 2); k++) {
  //     printf("%f ", ssums[i][k]);
  //   }
  // }
  // printf("\n");
  //
  // printf("snss: ");
  // for (int i = 0; i < 3; i++) {
  //   printf("%f ", snss[i]);
  // }
  // printf("\n");
  //
  // printf("schis: ");
  // for (int i = 0; i < 3; i++) {
  //   printf("%f ", schis[i]);
  // }
  // printf("\n");

  int return_code = 0;

  int j, j2, p, q, q2, q3, count, count2, start, end, mark, one = 1;
  double value, value2, sum, sum2, sumsq, mean, mean2, var, alpha, beta;

  int total, cflag, rflag;
  double scale, gc, sumhers, relax, likenull, like, like2, likeold, diff;
  double *thetas, *thetadiffs, *exps, *exps2, *jacks;
  double *sW, *sX, *sY, *sXTX, *sXTX2, *sXTY, *sXTXs, *sXTYs, *sT, *sTb, *sT2,
      *AI, *AI2, *AI3, *BI, *J, *JAI, *JAIJT;

  // assuming statj = c * (1 + nja + njvj b), where c is gif, a is intercept/n,
  // b is h2SNP/ssums easier to write as (statj-1) = nj/nvj cbn + (c-1) + nj/n
  // can = sT theta, where n is average nj

  // set total and maybe num_blocks
  total = num_parts + gcon + cept;
  if (num_blocks == -1 || num_blocks > length) {
    num_blocks = length;
  }

  // allocate variables

  thetas = malloc(sizeof(double) * total);
  thetadiffs = malloc(sizeof(double) * total);
  exps = malloc(sizeof(double) * length);
  exps2 = malloc(sizeof(double) * length);
  if (num_blocks != -9999) {
    jacks = malloc(sizeof(double) * (total + 1 + num_parts) * num_blocks);
  }

  sW = malloc(sizeof(double) * length);
  sX = malloc(sizeof(double) * total * length);
  sY = malloc(sizeof(double) * length);
  sXTX = malloc(sizeof(double) * total * total);
  sXTX2 = malloc(sizeof(double) * total * total);
  sXTY = malloc(sizeof(double) * total);
  sXTXs = malloc(sizeof(double) * total * total);
  sXTYs = malloc(sizeof(double) * total);

  sT = malloc(sizeof(double) * length * total);
  sTb = malloc(sizeof(double) * length * total);
  sT2 = malloc(sizeof(double) * length * total);
  AI = malloc(sizeof(double) * total * total);
  AI2 = malloc(sizeof(double) * total);
  AI3 = malloc(sizeof(double) * total * total);
  BI = malloc(sizeof(double) * total);
  J = malloc(sizeof(double) * total * total);
  JAI = malloc(sizeof(double) * total * total);
  JAIJT = malloc(sizeof(double) * total * total);

  // set variables that stay the same throughout

  // get scale - weighted average sample size
  sum = 0;
  sum2 = 0;
  for (j = 0; j < length; j++) {
    sum += snss[j] / stags[j];
    sum2 += pow(stags[j], -1);
  }
  for (j = 0; j < ncv; j++) {
    j2 = cvindex[j];
    sum -= snss[j2] / stags[j2];
    sum2 -= pow(stags[j2], -1);
  }
  scale = sum / sum2;

  // set sT and sTb (neither used if chisol=0 and use jackknifing for SDs)
  for (j = 0; j < length; j++) {
    for (q = 0; q < num_parts; q++) {
      sT[j + q * length] = snss[j] / scale * svars[q][j];
    }
    if (gcon == 1) {
      sT[j + num_parts * length] = 1;
    }
    if (cept == 1) {
      sT[j + (num_parts + gcon) * length] = snss[j] / scale;
    }
  }

  // sTb matches sT, except blanked out for cv predictors
  for (j = 0; j < length; j++) {
    for (q = 0; q < total; q++) {
      sTb[j + q * length] = sT[j + q * length];
    }
  }
  for (j = 0; j < ncv; j++) {
    for (q = 0; q < total; q++) {
      sTb[cvindex[j] + q * length] = 0;
    }
  }

  ////////

  // get null likelihood
  if (cept + gcon > 0) // null model will be E[Sj]=constant
  {
    sum = 0;
    sum2 = 0;
    for (j = 0; j < length; j++) {
      sum += schis[j] / stags[j];
      sum2 += pow(stags[j], -1);
    }
    for (j = 0; j < ncv; j++) {
      j2 = cvindex[j];
      sum -= schis[j2] / stags[j2];
      sum2 -= pow(stags[j2], -1);
    }
    value2 = sum / sum2;
  } else {
    value2 = 1;
  }

  sum = 0;
  sum2 = 0;
  value = 0;
  for (j = 0; j < length; j++) {
    sum += schis[j] / value2 / stags[j];
    sum2 += pow(stags[j], -1);
    value += log(schis[j] * value2) / stags[j];
  }
  for (j = 0; j < ncv; j++) {
    j2 = cvindex[j];
    sum -= schis[j2] / value2 / stags[j2];
    sum2 -= pow(stags[j2], -1);
    value -= log(schis[j2] * value2) / stags[j2];
  }
  likenull = -.5 * sum - .5 * value - .5 * sum2 * log(2 * M_PI);

  // set starting model

  if (sflag == 0 || sflag == 1) // starting model assumes no causal variation
                                // (plus gcon=1 and cept=1)
  {
    for (q = 0; q < total; q++) {
      thetas[q] = 0;
    }
  }

  if (sflag == 2 || sflag == 3 || sflag == 5) // starting model saved in stats
  {
    gc = 1;
    if (gcon == 1) {
      gc = stats[num_parts];
    }
    for (q = 0; q < num_parts; q++) {
      thetas[q] = stats[q] * gc / ssums[q][q] * scale;
    }
    if (gcon == 1) {
      thetas[num_parts] = (gc - 1);
    }
    if (cept == 1) {
      thetas[num_parts + gcon] = stats[num_parts + gcon] * gc;
    }
  }

  if (sflag == 4) // starting model copies ldsc - theta = (sum_j stat -1) /
                  // sum_jq (n_j var_jq)
  {
    value = 0;
    value2 = 0;
    for (j = 0; j < length; j++) {
      value += schis[j] - 1;
      for (q = 0; q < num_parts; q++) {
        value2 += snss[j] * svars[q][j];
      }
    }
    for (j = 0; j < ncv; j++) {
      value -= schis[cvindex[j]] - 1;
      for (q = 0; q < num_parts; q++) {
        value2 -= snss[cvindex[j]] * svars[q][cvindex[j]];
      }
    }

    for (q = 0; q < num_parts; q++) {
      thetas[q] = value / value2;
    }
    if (gcon == 1) {
      thetas[num_parts] = 0;
    }
    if (cept == 1) {
      thetas[num_parts + gcon] = 0;
    }
  }

  ////////

  // now iterate - rflag indicates type of move
  count = 0;
  cflag = 1;
  rflag = 1; // 0 if using multiNR (or least-squares), 1 if about to do
             // singleNR, 2 if just done singleNR
  while (1) {
    if (chisol == 0 || count == 0) // get exps and likelihood (already have for
                                   // NR solver after first iteration)
    {
      for (j = 0; j < length; j++) {
        exps[j] = 1;
      }
      alpha = 1.0;
      beta = 1.0;
      dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta,
             exps, &one);
      for (j = 0; j < length; j++) {
        if (exps[j] <= 0) {
          exps[j] = 1e-6;
        }
      }

      sum = 0;
      sum2 = 0;
      value = 0;
      for (j = 0; j < length; j++) {
        sum += schis[j] / exps[j] / stags[j];
        sum2 += pow(stags[j], -1);
        value += log(schis[j] * exps[j]) / stags[j];
      }
      for (j = 0; j < ncv; j++) {
        j2 = cvindex[j];
        sum -= schis[j2] / exps[j2] / stags[j2];
        sum2 -= pow(stags[j2], -1);
        value -= log(schis[j2] * exps[j2]) / stags[j2];
      }
      like = -.5 * sum - .5 * value - .5 * sum2 * log(2 * M_PI);
    }

    if (count > 0) // set diff
    {
      diff = like - likeold;
    }
    likeold = like;

    // print update
    gc = 1;
    if (gcon == 1) {
      gc = thetas[num_parts] + 1;
    }
    sumhers = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers += thetas[q] / gc * ssums[q][q] / scale;
    }

    // if (count == 0) {
    //   printf("Start\t");
    // } else {
    //   printf("%d\t", count);
    // }
    // printf("%.4f\t", sumhers);
    // if (gcon == 1) {
    //   printf("%.4f\t", gc);
    // }
    // if (cept == 1) {
    //   printf("%.4f\t", 1 + thetas[num_parts + gcon] / gc);
    // }
    // printf("%.2f\t", like);
    // if (count == 0) {
    //   printf("n/a\t\t%.6f\n", tol);
    // } else {
    //   printf("%.6f\t%.6f\n", diff, tol);
    // }

    if (sflag == 3) {
      break;
    } // only wanted expectations and likelihood
    if (sflag == 4 && num_parts > 1 && count == 1) {
      break;
    } // multi-tagging ldsc uses a single iteration
    if (count > 0) {
      if (fabs(diff) < tol && (rflag == 0 || chisol == 0)) {
        break;
      }
    }
    if (count == maxiter) {
      printf(
          "\nWarning, the optimizer failed to converge within %d iterations\n",
          maxiter);
      cflag = 0;
      break;
    }

    ////////

    // update thetas

    if (chisol == 0) // use least squares (will always accept proposed move)
    {
      // load up sW, sX and sY - regressing Sj-1 on nj*vjk, 1, nj
      for (j = 0; j < length; j++) {
        sW[j] = stags[j] * pow(exps[j], 2);
        value = pow(sW[j], -0.5);
        for (q = 0; q < num_parts; q++) {
          sX[j + q * length] = snss[j] / scale * svars[q][j] * value;
        }
        if (gcon == 1) {
          sX[j + num_parts * length] = value;
        }
        if (cept == 1) {
          sX[j + (num_parts + gcon) * length] = snss[j] / scale * value;
        }
        sY[j] = (schis[j] - 1) * value;
      }

      for (j = 0; j < ncv; j++) // blank out sX for cv predictors
      {
        for (q = 0; q < total; q++) {
          sX[cvindex[j] + q * length] = 0;
        }
      }

      // new theta is (sXTX)^-1 sXTY
      alpha = 1.0;
      beta = 0.0;
      dgemm_("T", "N", &total, &total, &length, &alpha, sX, &length, sX,
             &length, &beta, sXTX, &total);
      dgemv_("T", &length, &total, &alpha, sX, &length, sY, &one, &beta, sXTY,
             &one);

      for (q = 0; q < total; q++) {
        thetas[q] = sXTY[q];
      }
      int code = eigen_invert(sXTX, total, sXTX2, 1, thetas);
      if (code != 0) {
        printf("Error in eigen_invert\n");
        return_code = 1;
        goto cleanup;
      }
    } else // use either single nr or multi nr (only move if good for
           // likelihood)
    {
      if (rflag == 1) // single nr - for first iteration or if multi fails
      {
        for (q = 0; q < total; q++) {
          // get derivs for q - use sTb so as not to include contributions of cv
          // predictors
          value = 0;
          value2 = 0;
          for (j = 0; j < length; j++) {
            value += .5 * (schis[j] - exps[j]) / stags[j] * pow(exps[j], -2) *
                     sTb[j + q * length];
            value2 += (schis[j] - .5 * exps[j]) / stags[j] * pow(exps[j], -3) *
                      pow(sTb[j + q * length], 2);
          }

          // get proposed move
          thetadiffs[q] = value / value2;

          // move and get expectations
          thetas[q] += thetadiffs[q];

          for (j = 0; j < length; j++) {
            exps2[j] = 1;
          }
          alpha = 1.0;
          beta = 1.0;
          dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta,
                 exps2, &one);

          for (j = 0; j < length; j++) {
            if (exps2[j] <= 0) {
              exps2[j] = 1e-6;
            }
          }

          // get corresponding likelihood
          sum = 0;
          sum2 = 0;
          value = 0;
          for (j = 0; j < length; j++) {
            sum += schis[j] / exps2[j] / stags[j];
            sum2 += pow(stags[j], -1);
            value += log(schis[j] * exps2[j]) / stags[j];
          }
          for (j = 0; j < ncv; j++) {
            j2 = cvindex[j];
            sum -= schis[j2] / exps2[j2] / stags[j2];
            sum2 -= pow(stags[j2], -1);
            value -= log(schis[j2] * exps2[j2]) / stags[j2];
          }
          like2 = -.5 * sum - .5 * value - .5 * sum2 * log(2 * M_PI);

          if (like2 > like - tol) // accept move
          {
            like = like2;
            for (j = 0; j < length; j++) {
              exps[j] = exps2[j];
            }
          } else // move back
          {
            thetas[q] -= thetadiffs[q];
          }
        } // end of q loop
        rflag = 2;
      } else // multi nr - use BI and AI
      {
        // load up sT2=(stat-exps/2)*w/exps^3*sTb - will be blanked for
        // cvpredictors
        for (j = 0; j < length; j++) {
          value = (schis[j] - .5 * exps[j]) / stags[j] * pow(exps[j], -3);
          for (q = 0; q < total; q++) {
            sT2[j + q * length] = value * sTb[j + q * length];
          }
        }

        // get AI (-2nd deriv) and BI (1st deriv) - cvpredictors already taken
        // care of (by blanking sTb and ST2)
        alpha = 1.0;
        beta = 0.0;
        dgemm_("T", "N", &total, &total, &length, &alpha, sTb, &length, sT2,
               &length, &beta, AI, &total);
        for (q = 0; q < total; q++) {
          BI[q] = 0;
          for (j = 0; j < length; j++) {
            BI[q] += .5 * (schis[j] - exps[j]) / stags[j] * pow(exps[j], -2) *
                     sTb[j + q * length];
          }
        }

        int code = eigen_invert(AI, total, AI2, -1, AI3);
        if (code != 0) {
          printf("Error in eigen_invert\n");
          return_code = 1;
          goto cleanup;
        }

        // get proposed move
        alpha = 1.0;
        beta = 0.0;
        dgemv_("N", &total, &total, &alpha, AI, &total, BI, &one, &beta,
               thetadiffs, &one);

        rflag = 1;
        relax = 1;
        while (relax > 0.0001) {
          // move relax*thetadiffs and get expectations
          for (q = 0; q < total; q++) {
            thetas[q] += relax * thetadiffs[q];
          }

          for (j = 0; j < length; j++) {
            exps2[j] = 1;
          }
          alpha = 1.0;
          beta = 1.0;
          dgemv_("N", &length, &total, &alpha, sT, &length, thetas, &one, &beta,
                 exps2, &one);

          for (j = 0; j < length; j++) {
            if (exps2[j] <= 0) {
              exps2[j] = 1e-6;
            }
          }

          // get corresponding likelihood
          sum = 0;
          sum2 = 0;
          value = 0;
          for (j = 0; j < length; j++) {
            sum += schis[j] / exps2[j] / stags[j];
            sum2 += pow(stags[j], -1);
            value += log(schis[j] * exps2[j]) / stags[j];
          }
          for (j = 0; j < ncv; j++) {
            j2 = cvindex[j];
            sum -= schis[j2] / exps2[j2] / stags[j2];
            sum2 -= pow(stags[j2], -1);
            value -= log(schis[j2] * exps2[j2]) / stags[j2];
          }
          like2 = -.5 * sum - .5 * value - .5 * sum2 * log(2 * M_PI);

          if (like2 > like - tol) // accept move
          {
            like = like2;
            for (j = 0; j < length; j++) {
              exps[j] = exps2[j];
            }
            rflag = 0;
            break;
          } else // move back and next turn try smaller move
          {
            for (q = 0; q < total; q++) {
              thetas[q] -= relax * thetadiffs[q];
            }
            relax *= .5;
          }
        }
      }
    } // end of NR

    count++;
  } // end of while loop

  // load up first column of stats which contain Q hers, gc, na, sum of Q hers,
  // Q cats
  gc = 1;
  if (gcon == 1) {
    gc = thetas[num_parts] + 1;
  }
  sumhers = 0;
  for (q = 0; q < num_parts; q++) {
    sumhers += thetas[q] / gc * ssums[q][q] / scale;
  }

  for (q = 0; q < num_parts; q++) {
    stats[q] = thetas[q] / gc * ssums[q][q] / scale;
  }
  if (gcon == 1) {
    stats[num_parts] = gc;
  }
  if (cept == 1) {
    stats[num_parts + gcon] = thetas[num_parts + gcon] / gc;
  }
  stats[total] = sumhers;
  for (q = 0; q < num_parts; q++) {
    stats[total + 1 + q] = 0;
    for (q2 = 0; q2 < num_parts; q2++) {
      stats[total + 1 + q] += thetas[q2] / gc * ssums[q][q2] / scale;
    }
  }

  ////////

  if (likes !=
      NULL) // save likelihoods - currently redundant if using cv predictors
  {
    // already have chisq likelihoods
    likes[0] = likenull;
    likes[1] = like;

    // normal likelihood using starting weights
    if (cept + gcon > 0) // for null fit E[Sj]=constant
    {
      sum = 0;
      sum2 = 0;
      for (j = 0; j < length; j++) {
        sum += schis[j] / stags[j];
        sum2 += pow(stags[j], -1);
      }
      for (j = 0; j < ncv; j++) {
        j2 = cvindex[j];
        sum -= schis[j2] / stags[j2];
        sum2 -= pow(stags[j2], -1);
      }
      value2 = sum / sum2;
    } else {
      value2 = 1;
    }

    sum = 0;
    sum2 = 0;
    for (j = 0; j < length; j++) {
      sum += pow(schis[j] - value2, 2) / stags[j];
      sum2 += pow(stags[j], -1);
    }
    for (j = 0; j < ncv; j++) {
      j2 = cvindex[j];
      sum -= pow(schis[j2] - value2, 2) / stags[j2];
      sum2 -= pow(stags[j2], -1);
    }
    likes[2] = -.5 * sum2 * (1 + log(2 * M_PI * sum / sum2));

    sum = 0;
    sum2 = 0;
    for (j = 0; j < length; j++) {
      sum += pow(schis[j] - exps[j], 2) / stags[j];
      sum2 += pow(stags[j], -1);
    }
    for (j = 0; j < ncv; j++) {
      j2 = cvindex[j];
      sum -= pow(schis[j2] - exps[j2], 2) / stags[j2];
      sum2 -= pow(stags[j2], -1);
    }
    likes[3] = -.5 * sum2 * (1 + log(2 * M_PI * sum / sum2));

    // normal likelihood using final weights
    for (j = 0; j < length; j++) {
      sW[j] = stags[j] * pow(exps[j], 2);
    }

    if (cept + gcon > 0) // for null fit E[Sj]=constant
    {
      sum = 0;
      sum2 = 0;
      for (j = 0; j < length; j++) {
        sum += schis[j] / sW[j];
        sum2 += pow(sW[j], -1);
      }
      for (j = 0; j < ncv; j++) {
        j2 = cvindex[j];
        sum -= schis[j2] / sW[j2];
        sum2 -= pow(sW[j2], -1);
      }
      value2 = sum / sum2;
    } else {
      value2 = 1;
    }

    sum = 0;
    sum2 = 0;
    for (j = 0; j < length; j++) {
      sum += pow(schis[j] - value2, 2) / sW[j];
      sum2 += pow(sW[j], -1);
    }
    for (j = 0; j < ncv; j++) {
      j2 = cvindex[j];
      sum -= pow(schis[j2] - value2, 2) / sW[j2];
      sum2 -= pow(sW[j2], -1);
    }
    likes[4] = -.5 * sum2 * (1 + log(2 * M_PI * sum / sum2));

    sum = 0;
    sum2 = 0;
    for (j = 0; j < length; j++) {
      sum += pow(schis[j] - exps[j], 2) / sW[j];
      sum2 += pow(sW[j], -1);
    }
    for (j = 0; j < ncv; j++) {
      j2 = cvindex[j];
      sum -= pow(schis[j2] - exps[j2], 2) / sW[j2];
      sum2 -= pow(sW[j2], -1);
    }
    likes[5] = -.5 * sum2 * (1 + log(2 * M_PI * sum / sum2));

    // save cflag
    likes[6] = cflag;
  }

  // get negative counts
  count = 0;
  count2 = 0;
  sum = 0;
  sum2 = 0;
  for (j = 0; j < length; j++) {
    value = 0;
    for (q = 0; q < num_parts; q++) {
      value += svars[q][j] * thetas[q];
    }
    if (value < 0) {
      count++;
      sum += pow(stags[j], -1);
    }
    value2 = 1 + snss[j] / scale * value;
    if (gcon == 1) {
      value2 += thetas[num_parts];
    }
    if (cept == 1) {
      value2 += snss[j] * thetas[num_parts + gcon] / scale;
    }
    if (value2 <= 0) {
      count2++;
      sum2 += pow(stags[j], -1);
    }
  }

  if (likes != NULL) // save them
  {
    likes[7] = count;
    likes[8] = sum;
    likes[9] = count2;
    likes[10] = sum2;
  }

  if (ncv > 0) // save expectations for cv predictors
  {
    for (j = 0; j < ncv; j++) {
      cvexps[j] = exps[cvindex[j]];
    }
  }

  ////////

  if (ncv == 0 &&
      (sflag == 0 || sflag == 2 || sflag == 4 || sflag == 5)) // get SDs
  {
    if (num_blocks == -9999) // get from second derivative of likelihood (still
                             // valid to use thetas)
    {
      // get inverse AI for current state
      for (j = 0; j < length; j++) {
        value = (schis[j] - .5 * exps[j]) / stags[j] * pow(exps[j], -3);
        for (q = 0; q < total; q++) {
          sT2[j + q * length] = value * sTb[j + q * length];
        }
      }

      alpha = 1.0;
      beta = 0.0;
      dgemm_("T", "N", &total, &total, &length, &alpha, sTb, &length, sT2,
             &length, &beta, AI, &total);
      int code = eigen_invert(AI, total, AI2, -1, AI3);
      if (code != 0) {
        printf("Error in eigen_invert\n");
        return_code = 1;
        goto cleanup;
      }

      // AI provides variances for thetas = (cb*scale, (c-1), ca*scale) - for
      // transformed variances, must compute J invAI JT, where Jij=dnewi/dthetaj

      // make sure we have gc
      gc = 1;
      if (gcon == 1) {
        gc = thetas[num_parts] + 1;
      }

      // first want variances of (hers, c, na) = (b ssums, c, a*scale) =
      // (t1*ssums/(t2+1)/scale, t2+1, t3/(t2+1))
      for (q = 0; q < total; q++) {
        for (q2 = 0; q2 < total; q2++) {
          J[q + q2 * total] = 0;
        }
      }
      for (q = 0; q < num_parts; q++) {
        J[q + q * total] = ssums[q][q] / gc / scale;
      }
      if (gcon == 1) {
        for (q = 0; q < num_parts; q++) {
          J[q + num_parts * total] =
              -thetas[q] * ssums[q][q] / scale * pow(gc, -2);
        }
        J[num_parts + num_parts * total] = 1;
      }
      if (cept == 1) {
        if (gcon == 1) {
          J[num_parts + gcon + num_parts * total] =
              -thetas[num_parts + gcon] * pow(gc, -2);
        }
        J[num_parts + gcon + (num_parts + gcon) * total] = 1.0 / gc;
      }

      alpha = 1.0;
      beta = 0.0;
      dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total,
             &beta, JAI, &total);
      dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total,
             &beta, JAIJT, &total);

      for (q = 0; q < total; q++) {
        if (JAIJT[q + q * total] >= 0) {
          stats[q + total + 1 + num_parts] = pow(JAIJT[q + q * total], .5);
        } else {
          stats[q + total + 1 + num_parts] = -9999;
        }
      }

      if (cohers != NULL) {
        for (q = 0; q < num_parts; q++) {
          for (q2 = 0; q2 < num_parts; q2++) {
            cohers[q + q2 * num_parts] = JAIJT[q + q2 * total];
          }
        }
      }

      // can also get variance of sumhers
      sum = 0;
      for (q = 0; q < num_parts; q++) {
        for (q2 = 0; q2 < num_parts; q2++) {
          sum += JAIJT[q + q2 * total];
        }
      }
      if (sum >= 0) {
        stats[total + total + 1 + num_parts] = pow(sum, .5);
      } else {
        stats[total + total + 1 + num_parts] = -9999;
      }

      // and variances of category heritabilities
      for (q = 0; q < num_parts; q++) {
        sum = 0;
        for (q2 = 0; q2 < num_parts; q2++) {
          for (q3 = 0; q3 < num_parts; q3++) {
            sum += JAIJT[q2 + q3 * total] * ssums[q][q2] / ssums[q2][q2] *
                   ssums[q][q3] / ssums[q3][q3];
          }
        }
        if (sum >= 0) {
          stats[total + 1 + q + total + 1 + num_parts] = pow(sum, .5);
        } else {
          stats[total + 1 + q + total + 1 + num_parts] = -9999;
        }
      }

      // next want variances of (shares, c, ca) = (b ssums / sum b ssums,
      // c/avsum, ca) - last two irrelevant
      sum = 0;
      for (q = 0; q < num_parts; q++) {
        sum += thetas[q] * ssums[q][q];
      }
      for (q = 0; q < total; q++) {
        for (q2 = 0; q2 < total; q2++) {
          J[q + q2 * total] = 0;
        }
      }
      for (q = 0; q < num_parts; q++) {
        for (q2 = 0; q2 < num_parts; q2++) {
          J[q + q2 * total] =
              -thetas[q] * ssums[q][q] * ssums[q2][q2] * pow(sum, -2);
        }
        J[q + q * total] += ssums[q][q] / sum;
      }
      if (gcon == 1) {
        J[num_parts + num_parts * total] = 1;
      }
      if (cept == 1) {
        J[num_parts + gcon + (num_parts + gcon) * total] = 1;
      }

      alpha = 1.0;
      beta = 0.0;
      dgemm_("N", "N", &total, &total, &total, &alpha, J, &total, AI, &total,
             &beta, JAI, &total);
      dgemm_("N", "T", &total, &total, &total, &alpha, JAI, &total, J, &total,
             &beta, JAIJT, &total);

      for (q = 0; q < total; q++) {
        if (JAIJT[q + q * total] >= 0) {
          stats[q + 2 * (total + 1 + num_parts)] =
              pow(JAIJT[q + q * total], .5);
        } else {
          stats[q + 2 * (total + 1 + num_parts)] = -9999;
        }
      }

      // and variances of enrichments
      for (q = 0; q < num_parts; q++) {
        sum = 0;
        for (q2 = 0; q2 < num_parts; q2++) {
          for (q3 = 0; q3 < num_parts; q3++) {
            sum += JAIJT[q2 + q3 * total] * ssums[q][q2] / ssums[q2][q2] *
                   ssums[q][q3] / ssums[q3][q3];
          }
        }
        if (sum >= 0) {
          stats[total + 1 + q + 2 * (total + 1 + num_parts)] = pow(sum, .5);
        } else {
          stats[total + 1 + q + 2 * (total + 1 + num_parts)] = -9999;
        }
      }
    } else // jackknifing, using non-iterative least squares
    {
      // save sXTX and sXTY based on final expectations
      for (j = 0; j < length; j++) {
        sW[j] = stags[j] * pow(exps[j], 2);
        for (q = 0; q < num_parts; q++) {
          sX[j + q * length] = snss[j] / scale * svars[q][j] * pow(sW[j], -0.5);
        }
        if (gcon == 1) {
          sX[j + num_parts * length] = pow(sW[j], -0.5);
        }
        if (cept == 1) {
          sX[j + (num_parts + gcon) * length] =
              snss[j] / scale * pow(sW[j], -0.5);
        }
        sY[j] = (schis[j] - 1) * pow(sW[j], -0.5);
      }

      alpha = 1.0;
      beta = 0.0;
      dgemm_("T", "N", &total, &total, &length, &alpha, sX, &length, sX,
             &length, &beta, sXTX, &total);
      dgemv_("T", &length, &total, &alpha, sX, &length, sY, &one, &beta, sXTY,
             &one);

      for (q = 0; q < total; q++) {
        for (q2 = 0; q2 < total; q2++) {
          sXTXs[q + q2 * total] = sXTX[q + q2 * total];
        }
        sXTYs[q] = sXTY[q];
      }

      for (p = 0; p < num_blocks; p++) {
        start = (double)(p) / num_blocks * length;
        end = (double)(p + 1) / num_blocks * length;

        // reset sXTY and sXTY, then subtract and solve
        for (q = 0; q < total; q++) {
          for (q2 = 0; q2 < total; q2++) {
            sXTX[q + q2 * total] = sXTXs[q + q2 * total];
          }
          sXTY[q] = sXTYs[q];
        }

        alpha = -1.0;
        beta = 1.0;
        count = end - start;
        dgemm_("T", "N", &total, &total, &count, &alpha, sX + start, &length,
               sX + start, &length, &beta, sXTX, &total);
        dgemv_("T", &count, &total, &alpha, sX + start, &length, sY + start,
               &one, &beta, sXTY, &one);
        for (q = 0; q < total; q++) {
          thetas[q] = sXTY[q];
        }
        int code = eigen_invert(sXTX, total, sXTX2, 1, thetas);
        if (code != 0) {
          printf("Error in eigen_invert\n");
          return_code = 1;
          goto cleanup;
        }

        // save
        gc = 1;
        if (gcon == 1) {
          gc = thetas[num_parts] + 1;
        }
        sumhers = 0;
        for (q = 0; q < num_parts; q++) {
          sumhers += thetas[q] / gc * ssums[q][q] / scale;
        }

        mark = p * (total + 1 + num_parts);
        for (q = 0; q < num_parts; q++) {
          jacks[q + mark] = thetas[q] / gc * ssums[q][q] / scale;
        }
        if (gcon == 1) {
          jacks[num_parts + mark] = gc;
        }
        if (cept == 1) {
          jacks[num_parts + gcon + mark] = thetas[num_parts + gcon] / gc;
        }
        jacks[total + mark] = sumhers;
        for (q = 0; q < num_parts; q++) {
          jacks[total + 1 + q + mark] = 0;
          for (q2 = 0; q2 < num_parts; q2++) {
            jacks[total + 1 + q + mark] +=
                thetas[q2] / gc * ssums[q][q2] / scale;
          }
        }
      } // end of p loop

      // get sds for all stats (as is)
      for (q = 0; q < total + 1 + num_parts; q++) {
        sum = 0;
        sumsq = 0;
        for (p = 0; p < num_blocks; p++) {
          mark = p * (total + 1 + num_parts);
          sum += jacks[q + mark];
          sumsq += pow(jacks[q + mark], 2);
        }
        mean = sum / num_blocks;
        var = (num_blocks - 1) * (sumsq / num_blocks - pow(mean, 2));
        stats[q + total + 1 + num_parts] = pow(var, .5);
      }

      // now sds after dividing by sum of hers (for annotation shares,
      // considered instead dividing by base)
      for (q = 0; q < total + 1 + num_parts; q++) {
        sum = 0;
        sumsq = 0;
        for (p = 0; p < num_blocks; p++) {
          mark = p * (total + 1 + num_parts);
          sum += jacks[q + mark] / jacks[total + mark];
          sumsq += pow(jacks[q + mark] / jacks[total + mark], 2);
        }
        mean = sum / num_blocks;
        var = (num_blocks - 1) * (sumsq / num_blocks - pow(mean, 2));
        stats[q + 2 * (total + 1 + num_parts)] = pow(var, .5);
      }

      if (cohers !=
          NULL) // now coheritabilities - will be duplication of variances
      {
        for (q = 0; q < num_parts; q++) {
          for (q2 = 0; q2 < num_parts; q2++) {
            sum = 0;
            sum2 = 0;
            sumsq = 0;
            for (p = 0; p < num_blocks; p++) {
              mark = p * (total + 1 + num_parts);
              sum += jacks[q + mark];
              sum2 += jacks[q2 + mark];
              sumsq += jacks[q + mark] * jacks[q2 + mark];
            }
            mean = sum / num_blocks;
            mean2 = sum2 / num_blocks;
            var = (num_blocks - 1) * (sumsq / num_blocks - mean * mean2);
            cohers[q + q2 * num_parts] = var;
          }
        }
      }
    } // end of jackknifing
  } // end of getting SDs

  if (influs != NULL) // get factors for scaling hers to influences (can no
                      // longer use thetas)
  {
    // first get sum (stat/c-1-nja)^2/stags, weighted sum sq of test stats under
    // null (allowing for gc/cept)
    gc = 1;
    if (gcon == 1) {
      gc = stats[num_parts];
    }
    value = 0;
    if (cept == 1) {
      value = stats[num_parts + gcon] / scale;
    }
    sumsq = 0;
    for (j = 0; j < length; j++) {
      sumsq += pow(schis[j] / gc - 1 - snss[j] * value, 2) / stags[j];
    }

    for (q = 0; q < num_parts;
         q++) // get sum (stat/c-1-nja)(njvj)/stags /sumsq /ssums[q][q]
    {
      sum = 0;
      for (j = 0; j < length; j++) {
        sum += (schis[j] / gc - 1 - snss[j] * value) * snss[j] * svars[q][j] /
               stags[j];
      }
      influs[q] = sum / sumsq / ssums[q][q];
    }
  }

cleanup:
  free(thetas);
  free(thetadiffs);
  free(exps);
  free(exps2);
  if (num_blocks != -9999) {
    free(jacks);
  }
  free(sW);
  free(sX);
  free(sY);
  free(sXTX);
  free(sXTX2);
  free(sXTY);
  free(sXTXs);
  free(sXTYs);
  free(sT);
  free(sTb);
  free(sT2);
  free(AI);
  free(AI2);
  free(AI3);
  free(BI);
  free(J);
  free(JAI);
  free(JAIJT);
  return return_code;
}

int solve_cors(double *stats, int num_parts, int gcon, int cept, int num_blocks,
               int length, double *stags, double **svars, double **ssums,
               double *snss, double *schis, double *srhos, double *snss2,
               double *schis2, double *srhos2, double tol, int maxiter) {

  int return_code = 0;

  int j, q, q2, p, count, start, end, mark, one = 1;
  double sum, sum2, sum3, sumsq, mean, var, alpha, beta;

  int total, total2, total3, total4;
  double scale, scale2, scale3, gc, gc2, gc3, sumold, diff, sumhers, sumhers2,
      sumhers3, *snss3, *schis3, *exps, *exps2, *exps3;
  double *sW, *sX, *sY, *sXTX, *sXTX2, *sXTY, *sXTXs, *sXTYs, *thetas, *jacks;

  // set totals and maybe num_blocks
  total = 2 * (num_parts + gcon + cept) + num_parts + 1;
  total2 = num_parts + gcon + cept;
  total3 = num_parts + 1;
  total4 = num_parts + 2;
  if (num_blocks == -1 || num_blocks > length) {
    num_blocks = length;
  }

  // allocate variables

  snss3 = malloc(sizeof(double) * length);
  schis3 = malloc(sizeof(double) * length);
  exps = malloc(sizeof(double) * length);
  exps2 = malloc(sizeof(double) * length);
  exps3 = malloc(sizeof(double) * length);

  sW = malloc(sizeof(double) * length);
  sX = malloc(sizeof(double) * length * total4);
  sY = malloc(sizeof(double) * length);
  sXTX = malloc(sizeof(double) * total4 * total4);
  sXTX2 = malloc(sizeof(double) * total4);
  sXTY = malloc(sizeof(double) * total4);
  sXTXs = malloc(sizeof(double) * total4 * total4);
  sXTYs = malloc(sizeof(double) * total4);
  thetas = malloc(sizeof(double) * total4);
  jacks = malloc(sizeof(double) * (total + 4 + num_parts) * num_blocks);

  // set variables that stay the same throughout

  // snss3 = root (snss*snss2)
  for (j = 0; j < length; j++) {
    snss3[j] = pow(snss[j], .5) * pow(snss2[j], .5);
  }

  // product of z-statistics
  for (j = 0; j < length; j++) {
    if (srhos[j] * srhos2[j] > 0) {
      schis3[j] = pow(schis[j], .5) * pow(schis2[j], .5);
    } else {
      schis3[j] = -pow(schis[j], .5) * pow(schis2[j], .5);
    }
  }

  // get scale, scale2 and scale3 - weighted average sample sizes
  sum = 0;
  sum2 = 0;
  sum3 = 0;
  for (j = 0; j < length; j++) {
    sum += snss[j] / stags[j];
    sum2 += snss2[j] / stags[j];
    sum3 += pow(stags[j], -1);
  }
  scale = sum / sum3;
  scale2 = sum2 / sum3;
  scale3 = pow(scale, .5) * pow(scale2, .5);

  // null model blank
  sumold = 0;
  for (j = 0; j < length; j++) {
    exps[j] = 1;
  }

  count = 0;
  while (1) {
    // load up sW, sX and sY
    for (j = 0; j < length; j++) {
      sW[j] = stags[j] * pow(exps[j], 2);
      if (sW[j] <= 0) {
        sW[j] = 1e-6;
      }
      for (q = 0; q < num_parts; q++) {
        sX[j + q * length] = snss[j] / scale * svars[q][j] * pow(sW[j], -.5);
      }
      if (gcon == 1) {
        sX[j + num_parts * length] = pow(sW[j], -.5);
      }
      if (cept == 1) {
        sX[j + (num_parts + gcon) * length] = snss[j] / scale * pow(sW[j], -.5);
      }
      sY[j] = (schis[j] - 1) * pow(sW[j], -.5);
    }

    // get sXTX and sXTY
    alpha = 1.0;
    beta = 0.0;
    dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX,
           &length, &beta, sXTX, &total2);
    dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTY,
           &one);

    // solve for all parameters
    for (q = 0; q < total2; q++) {
      thetas[q] = sXTY[q];
    }
    int code = eigen_invert(sXTX, total2, sXTX2, 1, thetas);
    if (code != 0) {
      printf("Error in eigen_invert\n");
      return_code = 1;
      goto cleanup;
    }

    // get diff
    gc = 1;
    if (gcon == 1) {
      gc = thetas[num_parts] + 1;
    }
    sumhers = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers += thetas[q] / gc * ssums[q][q] / scale;
    }
    diff = sumhers - sumold;
    sumold = sumhers;

    // update exps - can use sX, but remember to multiply by sW[j] and add on
    // one
    alpha = 1.0;
    beta = 0.0;
    dgemv_("N", &length, &total2, &alpha, sX, &length, thetas, &one, &beta,
           exps, &one);
    for (j = 0; j < length; j++) {
      exps[j] = 1 + exps[j] * pow(sW[j], .5);
    }

    // see if breaking
    if (fabs(diff) < tol) {
      break;
    }
    if (count == maxiter) {
      break;
    }

    count++;
  }

  // put heritabilities, their sum, gc and cept into first column of stats (sum
  // goes right at end)
  gc = 1;
  if (gcon == 1) {
    gc = thetas[num_parts] + 1;
  }
  sumhers = 0;
  for (q = 0; q < num_parts; q++) {
    sumhers += thetas[q] / gc * ssums[q][q] / scale;
  }
  for (q = 0; q < num_parts; q++) {
    stats[q] = thetas[q] / gc * ssums[q][q] / scale;
  }
  if (gcon == 1) {
    stats[num_parts] = gc;
  }
  if (cept == 1) {
    stats[num_parts + gcon] = thetas[num_parts + gcon] / gc;
  }
  stats[total] = sumhers;

  // jackknife
  alpha = 1.0;
  beta = 0.0;
  dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX, &length,
         &beta, sXTXs, &total2);
  dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTYs,
         &one);

  for (p = 0; p < num_blocks; p++) {
    start = (double)p / num_blocks * length;
    end = (double)(p + 1) / num_blocks * length;
    count = end - start;

    // reset sXTX and sXTY
    for (q = 0; q < total2; q++) {
      for (q2 = 0; q2 < total2; q2++) {
        sXTX[q + q2 * total2] = sXTXs[q + q2 * total2];
      }
      sXTY[q] = sXTYs[q];
    }

    // then subtract and solve
    alpha = -1.0;
    beta = 1.0;
    dgemm_("T", "N", &total2, &total2, &count, &alpha, sX + start, &length,
           sX + start, &length, &beta, sXTX, &total2);
    dgemv_("T", &count, &total2, &alpha, sX + start, &length, sY + start, &one,
           &beta, sXTY, &one);
    for (q = 0; q < total2; q++) {
      thetas[q] = sXTY[q];
    }
    int code = eigen_invert(sXTX, total2, sXTX2, 1, thetas);
    if (code != 0) {
      printf("Error in eigen_invert\n");
      return_code = 1;
      goto cleanup;
    }

    // put heritabilities, their sum, gc and cept into jacks
    mark = p * (total + 4 + num_parts);
    gc = 1;
    if (gcon == 1) {
      gc = thetas[num_parts] + 1;
    }
    sumhers = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers += thetas[q] / gc * ssums[q][q] / scale;
    }
    for (q = 0; q < num_parts; q++) {
      jacks[q + mark] = thetas[q] / gc * ssums[q][q] / scale;
    }
    if (gcon == 1) {
      jacks[num_parts + mark] = gc;
    }
    if (cept == 1) {
      jacks[num_parts + gcon + mark] = thetas[num_parts + gcon] / gc;
    }
    jacks[total + mark] = sumhers;
  }

  ////////

  // null model blank
  sumold = 0;
  for (j = 0; j < length; j++) {
    exps2[j] = 1;
  }

  count = 0;
  while (1) {
    // load up sW, sX and sY
    for (j = 0; j < length; j++) {
      sW[j] = stags[j] * pow(exps2[j], 2);
      if (sW[j] <= 0) {
        sW[j] = 1e-6;
      }
      for (q = 0; q < num_parts; q++) {
        sX[j + q * length] = snss2[j] / scale2 * svars[q][j] * pow(sW[j], -.5);
      }
      if (gcon == 1) {
        sX[j + num_parts * length] = pow(sW[j], -.5);
      }
      if (cept == 1) {
        sX[j + (num_parts + gcon) * length] =
            snss2[j] / scale2 * pow(sW[j], -.5);
      }
      sY[j] = (schis2[j] - 1) * pow(sW[j], -.5);
    }

    // get sXTX and sXTY
    alpha = 1.0;
    beta = 0.0;
    dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX,
           &length, &beta, sXTX, &total2);
    dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTY,
           &one);

    // solve for all parameters
    for (q = 0; q < total2; q++) {
      thetas[q] = sXTY[q];
    }
    int code = eigen_invert(sXTX, total2, sXTX2, 1, thetas);
    if (code != 0) {
      printf("Error in eigen_invert\n");
      return_code = 1;
      goto cleanup;
    }

    // get diff
    gc2 = 1;
    if (gcon == 1) {
      gc2 = thetas[num_parts] + 1;
    }
    sumhers2 = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers2 += thetas[q] / gc2 * ssums[q][q] / scale2;
    }
    diff = sumhers2 - sumold;
    sumold = sumhers2;

    // update exps2 - can use sX, but remember to multiply by sW[j] and add on
    // one
    alpha = 1.0;
    beta = 0.0;
    dgemv_("N", &length, &total2, &alpha, sX, &length, thetas, &one, &beta,
           exps2, &one);
    for (j = 0; j < length; j++) {
      exps2[j] = 1 + exps2[j] * pow(sW[j], .5);
    }

    // see if breaking
    if (fabs(diff) < tol) {
      break;
    }
    if (count == maxiter) {
      break;
    }

    count++;
  }

  // put heritabilities, their sum, gc and cept into first column of stats (sum
  // goes right at end)
  gc2 = 1;
  if (gcon == 1) {
    gc2 = thetas[num_parts] + 1;
  }
  sumhers2 = 0;
  for (q = 0; q < num_parts; q++) {
    sumhers2 += thetas[q] / gc2 * ssums[q][q] / scale2;
  }
  for (q = 0; q < num_parts; q++) {
    stats[total2 + q] = thetas[q] / gc2 * ssums[q][q] / scale2;
  }
  if (gcon == 1) {
    stats[total2 + num_parts] = gc2;
  }
  if (cept == 1) {
    stats[total2 + num_parts + gcon] = thetas[num_parts + gcon] / gc2;
  }
  stats[total + 1] = sumhers2;

  // jackknife
  alpha = 1.0;
  beta = 0.0;
  dgemm_("T", "N", &total2, &total2, &length, &alpha, sX, &length, sX, &length,
         &beta, sXTXs, &total2);
  dgemv_("T", &length, &total2, &alpha, sX, &length, sY, &one, &beta, sXTYs,
         &one);

  for (p = 0; p < num_blocks; p++) {
    start = (double)p / num_blocks * length;
    end = (double)(p + 1) / num_blocks * length;
    count = end - start;

    // reset sXTX and sXTY
    for (q = 0; q < total2; q++) {
      for (q2 = 0; q2 < total2; q2++) {
        sXTX[q + q2 * total2] = sXTXs[q + q2 * total2];
      }
      sXTY[q] = sXTYs[q];
    }

    // then subtract and solve
    alpha = -1.0;
    beta = 1.0;
    dgemm_("T", "N", &total2, &total2, &count, &alpha, sX + start, &length,
           sX + start, &length, &beta, sXTX, &total2);
    dgemv_("T", &count, &total2, &alpha, sX + start, &length, sY + start, &one,
           &beta, sXTY, &one);
    for (q = 0; q < total2; q++) {
      thetas[q] = sXTY[q];
    }
    int code = eigen_invert(sXTX, total2, sXTX2, 1, thetas);
    if (code != 0) {
      printf("Error in eigen_invert\n");
      return_code = 1;
      goto cleanup;
    }

    // put heritabilities, their sum, gc and cept into jacks
    mark = p * (total + 4 + num_parts);
    gc2 = 1;
    if (gcon == 1) {
      gc2 = thetas[num_parts] + 1;
    }
    sumhers2 = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers2 += thetas[q] / gc2 * ssums[q][q] / scale2;
    }
    for (q = 0; q < num_parts; q++) {
      jacks[total2 + q + mark] = thetas[q] / gc2 * ssums[q][q] / scale2;
    }
    if (gcon == 1) {
      jacks[total2 + num_parts + mark] = gc2;
    }
    if (cept == 1) {
      jacks[total2 + num_parts + gcon + mark] = thetas[num_parts + gcon] / gc2;
    }
    jacks[total + 1 + mark] = sumhers2;
  }

  ////////

  // null model blank
  sumold = 0;
  for (j = 0; j < length; j++) {
    exps3[j] = 0;
  }

  count = 0;
  while (1) {
    // load up sW, sX and sY
    for (j = 0; j < length; j++) {
      sW[j] = stags[j] * (exps[j] * exps2[j] + pow(exps3[j], 2));
      if (sW[j] <= 0) {
        sW[j] = 1e-6;
      }
      for (q = 0; q < num_parts; q++) {
        sX[j + q * length] = snss3[j] / scale3 * svars[q][j] * pow(sW[j], -.5);
      }
      sX[j + num_parts * length] = pow(sW[j], -.5);
      sY[j] = schis3[j] * pow(sW[j], -.5);
    }

    // get sXTX and sXTY
    alpha = 1.0;
    beta = 0.0;
    dgemm_("T", "N", &total3, &total3, &length, &alpha, sX, &length, sX,
           &length, &beta, sXTX, &total3);
    dgemv_("T", &length, &total3, &alpha, sX, &length, sY, &one, &beta, sXTY,
           &one);

    // solve for all parameters
    for (q = 0; q < total3; q++) {
      thetas[q] = sXTY[q];
    }
    int code = eigen_invert(sXTX, total3, sXTX2, 1, thetas);
    if (code != 0) {
      printf("Error in eigen_invert\n");
      return_code = 1;
      goto cleanup;
    }

    // get diff
    gc = 1;
    if (gcon == 1) {
      gc = stats[num_parts];
    }
    gc2 = 1;
    if (gcon == 1) {
      gc2 = stats[total2 + num_parts];
    }
    sumhers3 = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers3 += thetas[q] * pow(gc * gc2, -.5) * ssums[q][q] / scale3;
    }
    diff = sumhers3 - sumold;
    sumold = sumhers3;

    // update exps3 - can use sX, but remember to multiply by sW[j] (but no need
    // to add on one)
    alpha = 1.0;
    beta = 0.0;
    dgemv_("N", &length, &total3, &alpha, sX, &length, thetas, &one, &beta,
           exps3, &one);
    for (j = 0; j < length; j++) {
      exps3[j] = exps3[j] * pow(sW[j], .5);
    }

    // see if breaking
    if (fabs(diff) < tol) {
      break;
    }
    if (count == maxiter) {
      break;
    }

    count++;
  }

  // put coheritabilities, their sum and the weird term into first column of
  // stats (sum goes right at end)
  gc = 1;
  if (gcon == 1) {
    gc = stats[num_parts];
  }
  gc2 = 1;
  if (gcon == 1) {
    gc2 = stats[total2 + num_parts];
  }
  gc3 = pow(gc * gc2, .5);
  sumhers3 = 0;
  for (q = 0; q < num_parts; q++) {
    sumhers3 += thetas[q] / gc3 * ssums[q][q] / scale3;
  }
  for (q = 0; q < num_parts; q++) {
    stats[2 * total2 + q] = thetas[q] / gc3 * ssums[q][q] / scale3;
  }
  stats[2 * total2 + num_parts] = thetas[num_parts] / gc3;
  stats[total + 2] = sumhers3;

  // and now also correlations
  stats[total + 3] =
      stats[total + 2] * pow(stats[total] * stats[total + 1], -.5);
  for (q = 0; q < num_parts; q++) {
    stats[total + 4 + q] =
        stats[2 * total2 + q] * pow(stats[q] * stats[total2 + q], -.5);
  }

  // jackknife
  alpha = 1.0;
  beta = 0.0;
  dgemm_("T", "N", &total3, &total3, &length, &alpha, sX, &length, sX, &length,
         &beta, sXTXs, &total3);
  dgemv_("T", &length, &total3, &alpha, sX, &length, sY, &one, &beta, sXTYs,
         &one);

  for (p = 0; p < num_blocks; p++) {
    start = (double)p / num_blocks * length;
    end = (double)(p + 1) / num_blocks * length;
    count = end - start;

    // reset sXTX and sXTY
    for (q = 0; q < total3; q++) {
      for (q2 = 0; q2 < total3; q2++) {
        sXTX[q + q2 * total3] = sXTXs[q + q2 * total3];
      }
      sXTY[q] = sXTYs[q];
    }

    // then subtract and solve
    alpha = -1.0;
    beta = 1.0;
    dgemm_("T", "N", &total3, &total3, &count, &alpha, sX + start, &length,
           sX + start, &length, &beta, sXTX, &total3);
    dgemv_("T", &count, &total3, &alpha, sX + start, &length, sY + start, &one,
           &beta, sXTY, &one);
    for (q = 0; q < total3; q++) {
      thetas[q] = sXTY[q];
    }
    int code = eigen_invert(sXTX, total3, sXTX2, 1, thetas);
    if (code != 0) {
      printf("Error in eigen_invert\n");
      return_code = 1;
      goto cleanup;
    }

    // put coheritabilities, their sum and the weird term into jacks
    mark = p * (total + 4 + num_parts);
    gc = 1;
    if (gcon == 1) {
      gc = jacks[num_parts + mark];
    }
    gc2 = 1;
    if (gcon == 1) {
      gc2 = jacks[total2 + num_parts + mark];
    }
    gc3 = pow(gc * gc2, .5);
    sumhers3 = 0;
    for (q = 0; q < num_parts; q++) {
      sumhers3 += thetas[q] / gc3 * ssums[q][q] / scale3;
    }
    for (q = 0; q < num_parts; q++) {
      jacks[2 * total2 + q + mark] = thetas[q] / gc3 * ssums[q][q] / scale3;
    }
    jacks[2 * total2 + num_parts + mark] = thetas[num_parts] / gc3;
    jacks[total + 2 + mark] = sumhers3;

    // and now also correlations
    mark = p * (total + 4 + num_parts);
    jacks[total + 3 + mark] =
        jacks[total + 2 + mark] *
        pow(jacks[total + mark] * jacks[total + 1 + mark], -.5);
    for (q = 0; q < num_parts; q++) {
      jacks[total + 4 + q + mark] =
          jacks[2 * total2 + q + mark] *
          pow(jacks[q + mark] * jacks[total2 + q + mark], -.5);
    }
  }

  ////////

  // now get all sds
  for (q = 0; q < total + 4 + num_parts; q++) {
    sum = 0;
    sumsq = 0;
    for (p = 0; p < num_blocks; p++) {
      mark = p * (total + 4 + num_parts);
      sum += jacks[q + mark];
      sumsq += pow(jacks[q + mark], 2);
    }
    mean = sum / num_blocks;
    var = (num_blocks - 1) * (sumsq / num_blocks - pow(mean, 2));
    stats[q + total + 4 + num_parts] = pow(var, .5);
  }

  // printf("Trait 1 heritability: %.4f (%.4f)\n", stats[total],
  //        stats[total + total + 4 + num_parts]);
  // printf("Trait 2 heritability: %.4f (%.4f)\n", stats[total + 1],
  //        stats[total + 1 + total + 4 + num_parts]);
  // printf("Coheritability: %.4f (%.4f)\n", stats[total + 2],
  //        stats[total + 2 + total + 4 + num_parts]);
  // printf("Correlation: %.4f (%.4f)\n\n", stats[total + 3],
  //        stats[total + 3 + total + 4 + num_parts]);

cleanup:
  free(snss3);
  free(schis3);
  free(exps);
  free(exps2);
  free(exps3);
  free(sW);
  free(sX);
  free(sY);
  free(sXTX);
  free(sXTX2);
  free(sXTY);
  free(sXTXs);
  free(sXTYs);
  free(thetas);
  free(jacks);
  return return_code;
}
