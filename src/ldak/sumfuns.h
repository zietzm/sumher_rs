#pragma once

int solve_sums(double *stats, double *likes, double *cohers, double *influs,
               int num_parts, int gcon, int cept, int num_blocks, int length,
               int ncv, int *cvindex, double *cvexps, double *stags,
               double **svars, double **ssums, double *snss, double *schis,
               double tol, int maxiter, int chisol, int sflag, char *filename);

int solve_cors(double *stats, int num_parts, int gcon, int cept, int num_blocks,
               int length, double *stags, double **svars, double **ssums,
               double *snss, double *schis, double *srhos, double *snss2,
               double *schis2, double *srhos2, double tol, int maxiter,
               char *filename);
