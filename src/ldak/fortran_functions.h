#pragma once

extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc);
extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
                   int *lda, double *x, int *incx, double *beta, double *y,
                   int *incy);
extern void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern void dpotri_(char *uplo, int *n, double *a, int *lda, int *info);
extern void dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                    double *b, int *ldb, int *info);
extern void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv,
                    double *work, int *lwork, int *info);
extern void dsytri_(char *uplo, int *n, double *a, int *lda, int *ipiv,
                    double *work, int *info);
extern void dsytrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                    int *ipiv, double *b, int *ldb, int *info);
extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                   double *w, double *work, int *lwork, int *info);
