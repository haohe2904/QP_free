#include "mex.h"
#include "lapack.h"
#include <string.h>

/* 声明 Fortran dormqr，使用 mwSignedIndex */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgIdAndTxt("apply_q_mex:nrhs", "需要4个输入: A, qr, tau, trans");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("apply_q_mex:nlhs", "最多1个输出");

    /* 输入参数 */
    const mxArray *mxA    = prhs[0];   /* 矩阵 A (n×p) */
    const mxArray *mxQR   = prhs[1];   /* qr 矩阵 (n×n) */
    const mxArray *mxTAU  = prhs[2];   /* tau 向量 (n×1) */
    const mxArray *mxTRANS= prhs[3];   /* 转置选项 'N' 或 'T' */

    if (!mxIsDouble(mxA) || mxIsComplex(mxA) ||
        !mxIsDouble(mxQR) || mxIsComplex(mxQR) ||
        !mxIsDouble(mxTAU) || mxIsComplex(mxTAU))
        mexErrMsgIdAndTxt("apply_q_mex:type", "所有输入必须为双精度实数");

    /* 获取维度 */
    mwSignedIndex m = (mwSignedIndex)mxGetM(mxA);    /* 行数 = n */
    mwSignedIndex n = (mwSignedIndex)mxGetN(mxA);    /* 列数 = p */
    mwSignedIndex k = (mwSignedIndex)mxGetNumberOfElements(mxTAU); /* 反射器个数 = n */
    mwSignedIndex lda = (mwSignedIndex)mxGetM(mxQR); /* = n */
    mwSignedIndex ldc = m;

    /* 复制输入矩阵 A，因为 dormqr 会修改 C */
    double *C = (double*)mxMalloc(m * n * sizeof(double));
    if (!C) mexErrMsgIdAndTxt("apply_q_mex:mem", "内存分配失败");
    memcpy(C, mxGetPr(mxA), m * n * sizeof(double));

    /* 获取 QR 和 tau 数据指针 */
    double *QR = mxGetPr(mxQR);
    double *TAU = mxGetPr(mxTAU);

    /* 解析 trans 选项 */
    char trans;
    char trans_str[2];
    if (mxGetString(mxTRANS, trans_str, 2) != 0)
        mexErrMsgIdAndTxt("apply_q_mex:trans", "trans 必须是 'N' 或 'T'");
    trans = trans_str[0];
    if (trans != 'N' && trans != 'T')
        mexErrMsgIdAndTxt("apply_q_mex:trans", "trans 必须是 'N' 或 'T'");

    char side = 'L';   /* 左乘 Q */

    /* 查询最优工作区大小 */
    mwSignedIndex lwork = -1;
    double work_opt;
    mwSignedIndex info = 0;
    dormqr(&side, &trans, &m, &n, &k, QR, &lda, TAU, C, &ldc, &work_opt, &lwork, &info);
    if (info != 0) {
        mxFree(C);
        mexErrMsgIdAndTxt("apply_q_mex:query", "工作区查询失败, info = %d", (int)info);
    }

    /* 分配工作区 */
    lwork = (mwSignedIndex)work_opt;
    double *work = (double*)mxMalloc(lwork * sizeof(double));
    if (!work) {
        mxFree(C);
        mexErrMsgIdAndTxt("apply_q_mex:mem", "工作区内存分配失败");
    }

    /* 实际调用 dormqr */
    dormqr(&side, &trans, &m, &n, &k, QR, &lda, TAU, C, &ldc, work, &lwork, &info);
    if (info != 0) {
        mxFree(C);
        mxFree(work);
        mexErrMsgIdAndTxt("apply_q_mex:dormqr", "dormqr 执行失败, info = %d", (int)info);
    }

    /* 创建输出矩阵 */
    plhs[0] = mxCreateDoubleMatrix((mwSize)m, (mwSize)n, mxREAL);
    double *out = mxGetPr(plhs[0]);
    memcpy(out, C, m * n * sizeof(double));

    mxFree(C);
    mxFree(work);
}