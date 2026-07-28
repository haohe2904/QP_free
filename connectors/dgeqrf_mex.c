#include "mex.h"
#include "lapack.h"
#include <string.h>

/* 直接声明 Fortran dgeqrf，使用 mwSignedIndex */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 1) mexErrMsgIdAndTxt("dgeqrf_mex:nrhs","需要1个输入");
    if (nlhs > 2) mexErrMsgIdAndTxt("dgeqrf_mex:nlhs","最多2个输出");
    
    const mxArray *mxA = prhs[0];
    if (!mxIsDouble(mxA) || mxIsComplex(mxA))
        mexErrMsgIdAndTxt("dgeqrf_mex:type","A必须是双精度实数矩阵");
    
    mwSignedIndex m = (mwSignedIndex)mxGetM(mxA);
    mwSignedIndex n = (mwSignedIndex)mxGetN(mxA);
    mwSignedIndex lda = m;
    mwSignedIndex min_mn = (m < n) ? m : n;
    
    double *A = (double*)mxMalloc(m * n * sizeof(double));
    if (!A) mexErrMsgIdAndTxt("dgeqrf_mex:mem","内存分配失败");
    memcpy(A, mxGetPr(mxA), m * n * sizeof(double));
    
    double *tau = (double*)mxMalloc(min_mn * sizeof(double));
    if (!tau) { mxFree(A); mexErrMsgIdAndTxt("dgeqrf_mex:mem","tau内存失败"); }
    
    mwSignedIndex lwork = -1;
    double work_opt;
    mwSignedIndex info = 0;
    dgeqrf(&m, &n, A, &lda, tau, &work_opt, &lwork, &info);
    if (info != 0) { mxFree(A); mxFree(tau); mexErrMsgIdAndTxt("dgeqrf_mex:query","工作区查询失败"); }
    
    lwork = (mwSignedIndex)work_opt;
    double *work = (double*)mxMalloc(lwork * sizeof(double));
    if (!work) { mxFree(A); mxFree(tau); mexErrMsgIdAndTxt("dgeqrf_mex:mem","工作区失败"); }
    
    dgeqrf(&m, &n, A, &lda, tau, work, &lwork, &info);
    if (info != 0) { mxFree(A); mxFree(tau); mxFree(work); mexErrMsgIdAndTxt("dgeqrf_mex:lapack","dgeqrf执行失败"); }
    
    plhs[0] = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
    memcpy(mxGetPr(plhs[0]), A, m * n * sizeof(double));
    plhs[1] = mxCreateDoubleMatrix((mwSize)min_mn, 1, mxREAL);
    memcpy(mxGetPr(plhs[1]), tau, min_mn * sizeof(double));
    
    mxFree(A); mxFree(tau); mxFree(work);
}