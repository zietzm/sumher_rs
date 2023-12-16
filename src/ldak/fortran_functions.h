#pragma once

extern void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
                   double *alpha, double *a, int *lda, double *b, int *ldb,
                   double *beta, double *c, int *ldc);
extern void sgemm_();
extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *a,
                   int *lda, double *x, int *incx, double *beta, double *y,
                   int *incy);
extern void sgemv_();
extern void dsymm_();
extern void ssymm_();
extern void dsymv_();
extern void ssymv_();
extern void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern void spotrf_();
extern void dpotri_(char *uplo, int *n, double *a, int *lda, int *info);
extern void spotri_();
extern void dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                    double *b, int *ldb, int *info);
extern void spotrs_();
extern void dsytrf_(char *uplo, int *n, double *a, int *lda, int *ipiv,
                    double *work, int *lwork, int *info);
extern void ssytrf_();
extern void dsytri_(char *uplo, int *n, double *a, int *lda, int *ipiv,
                    double *work, int *info);
extern void ssytri_();
extern void dsytrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                    int *ipiv, double *b, int *ldb, int *info);
extern void ssytrs_();
extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                   double *w, double *work, int *lwork, int *info);
extern void ssyev_();
extern void dsyevx_();
extern void dgesvd_();
extern void dgesdd_();
extern void dtrmm_();
extern void strmm_();
extern void dtrmv_();
extern void strmv_();
extern void sspmv_();
extern void spptrf_();
extern void spptri_();
extern void spptrs_();
