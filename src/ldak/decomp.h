#pragma once

double cg_solve(double *mat, int length, double *mat2, double *mat3, int ncol,
                double tol);

double eigen_invert(double *mat, int length, double *mat2, int ncol,
                    double *mat3, int type);

void eigen_strip(double *mat, int start, int end, int total, double strip);

double ldlt_invert(double *mat, int length, int ncol, double *mat2, int *info,
                   int type);

double cholesky_invert(double *mat, int length, int ncol, double *mat2,
                       int *info, int type);
