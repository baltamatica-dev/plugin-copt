#include "coptmex.h"

extern int COPT_SearchParamAttr(copt_prob *prob, const char *name, int *p_type);

extern int COPT_GetPSDSolution(copt_prob* prob,
                               double* psdColValue,
                               double* psdRowSlack,
                               double* psdRowDual,
                               double* psdColDual);

extern int COPT_LoadConeProb(copt_prob* prob,
  int nCol,
  int nRow,
  int nFree,
  int nPositive,
  int nBox,
  int nCone,
  int nRotateCone,
  int nPrimalExp,
  int nDualExp,
  int nPrimalPow,
  int nDualPow,
  int nPSD,
  int nQObjElem,
  int iObjSense,
  double dObjConst,
  const double* colObj,
  const int* qObjRow,
  const int* qObjCol,
  const double* qObjElem,
  const int* colMatBeg,
  const int* colMatCnt,
  const int* colMatIdx,
  const double* colMatElem,
  const double* rowRhs,
  const double* boxLower,
  const double* boxUpper,
  const int* coneDim,
  const int* rotateConeDim,
  const int* primalPowDim,
  const int* dualPowDim,
  const double* primalPowAlpha,
  const double* dualPowAlpha,
  const int* psdDim,
  const char* colType,
  char const* const* colNames,
  char const* const* rowNames,
  char const* const* psdColNames,
  int* outRowMap);

/* Convert status code from integer to string */
static const char *COPTMEX_statusInt2Str(int status) {
  switch (status) {
  case 0:  // unstarted
    return COPTMEX_STATUS_UNSTARTED;
  case 1:  // optimal
    return COPTMEX_STATUS_OPTIMAL;
  case 2:  // infeasible
    return COPTMEX_STATUS_INFEASIBLE;
  case 3:  // unbounded
    return COPTMEX_STATUS_UNBOUNDED;
  case 4:  // inf_or_unb
    return COPTMEX_STATUS_INF_OF_UNB;
  case 5:  // numerical
    return COPTMEX_STATUS_NUMERICAL;
  case 6:  // nodelimit
    return COPTMEX_STATUS_NODELIMIT;
  case 8:  // timeout
    return COPTMEX_STATUS_TIMEOUT;
  case 9:  // unfinished
    return COPTMEX_STATUS_UNFINISHED;
  case 10: // interrupted
    return COPTMEX_STATUS_INTERRUPTED;
  default: // unknown
    return "unknown";
  }
}

/* Convert objective sense from string to integer */
static int COPTMEX_objsenStr2Int(char *objsen) {
  if (mystrcmp(objsen, "max") == 0 || mystrcmp(objsen, "maximize") == 0) {
    return COPT_MAXIMIZE;
  }
  return COPT_MINIMIZE;
}

/* Convert objective sense from integer to string */
static const char *COPTMEX_objsenInt2Str(int objsen) {
  if (objsen == COPT_MAXIMIZE) {
    return "max";
  }
  return "min";
}

/* Extract string from MEX */
static int COPTMEX_getString(const bxArray *in_array, char **out_str) {
  int retcode = 0;

  int bufflen = bxGetNumberOfElements(in_array) + 1;
  char *buffer = (char *) bxCalloc(bufflen, sizeof(char));
  if (!buffer) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }
  
  bxAsCStr(in_array, buffer, bufflen);
  *out_str = buffer;

exit_cleanup:
  return retcode;
}

/* Free string allocated by MEX */
static void COPTMEX_freeString(char **in_str) {
  if (in_str != NULL && *in_str != NULL) {
    bxFree(*in_str);
  }
  return;
}

/* Display error message */
void COPTMEX_errorMsg(int errcode, const char *errinfo) {
  char *errid = NULL;
  char *errtxt = NULL;

  switch (errcode) {
  case COPTMEX_ERROR_BAD_TYPE:
    errid  = "coptmex:BadType";
    errtxt = "Invalid type of '%s'.";
    break;
  case COPTMEX_ERROR_BAD_NAME:
    errid  = "coptmex:BadName";
    errtxt = "Invalid name of '%s'.";
    break;
  case COPTMEX_ERROR_BAD_DATA:
    errid  = "coptmex:BadData";
    errtxt = "Invalid data of '%s'.";
    break;
  case COPTMEX_ERROR_BAD_NUM:
    errid  = "coptmex:BadNum";
    errtxt = "Invalid number of elemtents in '%s'.";
    break;
  case COPTMEX_ERROR_BAD_API:
    errid  = "coptmex:BadAPI";
    errtxt = "%s.";
    break;
  default:
    return;
  }

  mexErrMsgIdAndTxt(errid, errtxt, errinfo);
}

/* Convert CSC matrix to COO matrix */
void COPTMEX_csc2coo(bxArray *q, int *qMatRow, int *qMatCol, double *qMatElem) {
  int ncol = bxGetN(q);
  mwIndex *jc = bxGetJc(q);
  mwIndex *ir = bxGetIr(q);
  double *val = bxGetDoubles(q);

  for (int i = 0; i < ncol; ++i) {
    int nColElem = (int) jc[i];
    int nColLast = (int) jc[i + 1];
    for (; nColElem < nColLast; ++nColElem) {
      qMatRow[nColElem] = (int) ir[nColElem];
      qMatCol[nColElem] = i;
      qMatElem[nColElem] = val[nColElem];
    }
  }

  return;
}

/* Convert COO matrix to CSC matrix */
int COPTMEX_coo2csc(int nQElem, int *qMatRow, int *qMatCol, double *qMatElem, bxArray *q) {
  int ncol = bxGetN(q);
  mwIndex *jc = bxGetJc(q);
  mwIndex *ir = bxGetIr(q);
  double *val = bxGetDoubles(q);

  int *colMatCnt = (int *) bxCalloc(ncol, sizeof(int));
  if (!colMatCnt) {
    return COPT_RETCODE_MEMORY;
  }

  for (int i = 0; i < nQElem; ++i) {
    colMatCnt[qMatCol[i]]++;
  }

  jc[0] = 0;
  for (int i = 1; i <= ncol; ++i) {
    jc[i] = jc[i - 1] + colMatCnt[i - 1];
  }

  for (int i = 0; i < nQElem; ++i) {
    int iCol = qMatCol[i];
    int iElem = (int) jc[iCol];

    ir[iElem] = qMatRow[i];
    val[iElem] = qMatElem[i];

    jc[iCol]++;
  }

  for (int i = 0, last = 0; i <= ncol; ++i) {
    int tmp = jc[i];
    jc[i] = last;
    last = tmp;
  }

  bxFree(colMatCnt);
  return COPT_RETCODE_OK;
}

#include "coptinit.c"

/* Display banner */
int COPTMEX_dispBanner(void) {
  int retcode = 0;
  char msgbuf[COPT_BUFFSIZE];

  COPTMEX_CALL(COPT_GetBanner(msgbuf, COPT_BUFFSIZE));
  bxPrintf("%s", msgbuf);

exit_cleanup:
  return retcode;
}

/* Extract version of COPT */
int COPTMEX_getVersion(bxArray **out_version) {
  int retcode = COPT_RETCODE_OK;
  bxArray *retver = NULL;
  coptmex_mversion coptver;

  COPTMEX_initVersion(&coptver);

  coptver.major = bxCreateDoubleMatrix(1, 1, bxREAL);
  coptver.minor = bxCreateDoubleMatrix(1, 1, bxREAL);
  coptver.technical = bxCreateDoubleMatrix(1, 1, bxREAL);
  if (!coptver.major || !coptver.minor || !coptver.technical) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  *bxGetDoubles(coptver.major) = COPT_VERSION_MAJOR;
  *bxGetDoubles(coptver.minor) = COPT_VERSION_MINOR;
  *bxGetDoubles(coptver.technical) = COPT_VERSION_TECHNICAL;

  retver = bxCreateStructMatrix(1, 1, 0, NULL);
  if (!retver) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // 'major'
  bxAddField(retver, COPTMEX_VERSION_MAJOR);
  bxSetField(retver, 0, COPTMEX_VERSION_MAJOR, coptver.major);
  // 'minor'
  bxAddField(retver, COPTMEX_VERSION_MINOR);
  bxSetField(retver, 0, COPTMEX_VERSION_MINOR, coptver.minor);
  // 'technical'
  bxAddField(retver, COPTMEX_VERSION_TECHNICAL);
  bxSetField(retver, 0, COPTMEX_VERSION_TECHNICAL, coptver.technical);

  *out_version = retver;

exit_cleanup:
  if (retcode != COPT_RETCODE_OK) {
    *out_version = NULL;
  }

  return retcode;
}

/* Extract objective sense */
int COPTMEX_getObjsen(const bxArray *in_objsen, int *out_objsen) {
  int retcode = COPT_RETCODE_OK;
  char *objsen = NULL;

  COPTMEX_CALL(COPTMEX_getString(in_objsen, &objsen));
  *out_objsen = COPTMEX_objsenStr2Int(objsen);

exit_cleanup:
  COPTMEX_freeString(&objsen);
  return retcode;
}

/* Extract LP solution */
static int COPTMEX_getLpResult(copt_prob *prob, bxArray **out_lpresult) {
  int retcode = COPT_RETCODE_OK;
  bxArray *lpResult = NULL;
  coptmex_clpsol csol;
  coptmex_mlpsol msol;

  COPTMEX_initCLpSol(&csol);
  COPTMEX_initMLpSol(&msol);

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &csol.nRow));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &csol.nCol));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_PSDCOLS, &csol.nPSD));
  COPTMEX_CALL(COPT_GetIntAttr(prob, "PSDLens", &csol.nPSDLen));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_PSDCONSTRS, &csol.nPSDConstr));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_QCONSTRS, &csol.nQConstr));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_LPSTATUS, &csol.nStatus));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASBASIS, &csol.hasBasis));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASLPSOL, &csol.hasLpSol));

  msol.status      = bxCreateString(COPTMEX_statusInt2Str(csol.nStatus));
  msol.simplexiter = bxCreateDoubleMatrix(1, 1, bxREAL);
  msol.barrieriter = bxCreateDoubleMatrix(1, 1, bxREAL);
  msol.solvingtime = bxCreateDoubleMatrix(1, 1, bxREAL);
  if (!msol.status || !msol.simplexiter || !msol.barrieriter || !msol.solvingtime) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  if (csol.hasLpSol) {
    msol.objval  = bxCreateDoubleMatrix(1, 1, bxREAL);
    msol.value   = bxCreateDoubleMatrix(csol.nCol, 1, bxREAL);
    msol.redcost = bxCreateDoubleMatrix(csol.nCol, 1, bxREAL);
    msol.slack   = bxCreateDoubleMatrix(csol.nRow, 1, bxREAL);
    msol.dual    = bxCreateDoubleMatrix(csol.nRow, 1, bxREAL);
    if (!msol.objval || !msol.value || !msol.redcost || !msol.slack || !msol.dual) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    csol.colValue = bxGetDoubles(msol.value);
    csol.colDual  = bxGetDoubles(msol.redcost);
    csol.rowSlack = bxGetDoubles(msol.slack);
    csol.rowDual  = bxGetDoubles(msol.dual);

    if (csol.nQConstr > 0) {
      msol.qcslack = bxCreateDoubleMatrix(csol.nQConstr, 1, bxREAL);
      if (!msol.qcslack) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      csol.qRowSlack = bxGetDoubles(msol.qcslack);
    }

    if (csol.nPSD > 0) {
      msol.psdcolvalue = bxCreateDoubleMatrix(csol.nPSDLen, 1, bxREAL);
      msol.psdcoldual = bxCreateDoubleMatrix(csol.nPSDLen, 1, bxREAL);
      if (!msol.psdcolvalue || !msol.psdcoldual) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      csol.psdColValue = bxGetDoubles(msol.psdcolvalue);
      csol.psdColDual = bxGetDoubles(msol.psdcoldual);
    }

    if (csol.nPSDConstr > 0) {
      msol.psdrowslack = bxCreateDoubleMatrix(csol.nPSDConstr, 1, bxREAL);
      msol.psdrowdual = bxCreateDoubleMatrix(csol.nPSDConstr, 1, bxREAL);
      if (!msol.psdrowslack || !msol.psdrowdual) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      csol.psdRowSlack = bxGetDoubles(msol.psdrowslack);
      csol.psdRowDual = bxGetDoubles(msol.psdrowdual);
    }
  }

  if (csol.nStatus == COPT_LPSTATUS_INFEASIBLE ||
      csol.nStatus == COPT_LPSTATUS_UNBOUNDED) {
    int iReqFarkasRay = 0;
    COPTMEX_CALL(COPT_GetIntParam(prob, COPT_INTPARAM_REQFARKASRAY, &iReqFarkasRay));

    if (iReqFarkasRay) {
      if (csol.nStatus == COPT_LPSTATUS_INFEASIBLE) {
        msol.farkas = bxCreateDoubleMatrix(csol.nRow, 1, bxREAL);
        if (!msol.farkas) {
          retcode = COPT_RETCODE_MEMORY;
          goto exit_cleanup;
        }
        csol.dualFarkas = bxGetDoubles(msol.farkas);
      }

      if (csol.nStatus == COPT_LPSTATUS_UNBOUNDED) {
        msol.ray = bxCreateDoubleMatrix(csol.nCol, 1, bxREAL);
        if (!msol.ray) {
          retcode = COPT_RETCODE_MEMORY;
          goto exit_cleanup;
        }
        csol.primalRay = bxGetDoubles(msol.ray);
      }
    }
  }

  if (csol.hasBasis) {
    msol.varbasis    = bxCreateDoubleMatrix(csol.nCol, 1, bxREAL);
    msol.constrbasis = bxCreateDoubleMatrix(csol.nRow, 1, bxREAL);
    if (!msol.varbasis || !msol.constrbasis) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    csol.colBasis = (int *) bxCalloc(csol.nCol, sizeof(int));
    csol.rowBasis = (int *) bxCalloc(csol.nRow, sizeof(int));
    if (!csol.colBasis || !csol.rowBasis) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_SIMPLEXITER, &csol.nSimplexIter));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_BARRIERITER, &csol.nBarrierIter));
  COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_SOLVINGTIME, &csol.dSolvingTime));

  if (csol.hasLpSol) {
    COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_LPOBJVAL, &csol.dObjVal));
    COPTMEX_CALL(COPT_GetLpSolution(prob, csol.colValue, csol.rowSlack,
                 csol.rowDual, csol.colDual));

    if (csol.nQConstr > 0) {
      COPTMEX_CALL(COPT_GetQConstrInfo(prob, COPT_DBLINFO_SLACK, csol.nQConstr, NULL, csol.qRowSlack));
    }

    if (csol.nPSD > 0) {
      COPTMEX_CALL(COPT_GetPSDSolution(prob, csol.psdColValue, NULL, NULL, csol.psdColDual));
    }

    if (csol.nPSDConstr > 0) {
      COPTMEX_CALL(COPT_GetPSDSolution(prob, NULL, csol.psdRowSlack, csol.psdRowDual, NULL));
    }
  }

  if (csol.dualFarkas != NULL) {
    COPTMEX_CALL(COPT_GetRowInfo(prob, COPT_DBLINFO_DUALFARKAS, csol.nRow, NULL, csol.dualFarkas));
  }

  if (csol.primalRay != NULL) {
    COPTMEX_CALL(COPT_GetColInfo(prob, COPT_DBLINFO_PRIMALRAY, csol.nCol, NULL, csol.primalRay));
  }

  if (csol.hasBasis) {
    COPTMEX_CALL(COPT_GetBasis(prob, csol.colBasis, csol.rowBasis));
  }

  *bxGetDoubles(msol.simplexiter) = csol.nSimplexIter;
  *bxGetDoubles(msol.barrieriter) = csol.nBarrierIter;
  *bxGetDoubles(msol.solvingtime) = csol.dSolvingTime;

  if (csol.hasLpSol) {
    *bxGetDoubles(msol.objval) = csol.dObjVal;
  }

  if (csol.hasBasis) {
    double *colBasis_data = bxGetDoubles(msol.varbasis);
    double *rowBasis_data = bxGetDoubles(msol.constrbasis);
    for (int i = 0; i < csol.nCol; ++i) {
      colBasis_data[i] = csol.colBasis[i];
    }
    for (int i = 0; i < csol.nRow; ++i) {
      rowBasis_data[i] = csol.rowBasis[i];
    }

    bxFree(csol.colBasis);
    bxFree(csol.rowBasis);
  }

  lpResult = bxCreateStructMatrix(1, 1, 0, NULL);
  if (!lpResult) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // 'status'
  bxAddField(lpResult, COPTMEX_RESULT_STATUS);
  bxSetField(lpResult, 0, COPTMEX_RESULT_STATUS, msol.status);
  // 'simplexiter'
  bxAddField(lpResult, COPTMEX_RESULT_SIMITER);
  bxSetField(lpResult, 0, COPTMEX_RESULT_SIMITER, msol.simplexiter);
  // 'barrieriter'
  bxAddField(lpResult, COPTMEX_RESULT_BARITER);
  bxSetField(lpResult, 0, COPTMEX_RESULT_BARITER, msol.barrieriter);
  // 'solvingtime'
  bxAddField(lpResult, COPTMEX_RESULT_SOLVETIME);
  bxSetField(lpResult, 0, COPTMEX_RESULT_SOLVETIME, msol.solvingtime);

  if (csol.hasLpSol) {
    // 'objval'
    bxAddField(lpResult, COPTMEX_RESULT_OBJVAL);
    bxSetField(lpResult, 0, COPTMEX_RESULT_OBJVAL, msol.objval);
    // 'x'
    bxAddField(lpResult, COPTMEX_RESULT_VALUE);
    bxSetField(lpResult, 0, COPTMEX_RESULT_VALUE, msol.value);
    // 'rc'
    bxAddField(lpResult, COPTMEX_RESULT_REDCOST);
    bxSetField(lpResult, 0, COPTMEX_RESULT_REDCOST, msol.redcost);
    // 'slack;
    bxAddField(lpResult, COPTMEX_RESULT_SLACK);
    bxSetField(lpResult, 0, COPTMEX_RESULT_SLACK, msol.slack);
    // 'pi'
    bxAddField(lpResult, COPTMEX_RESULT_DUAL);
    bxSetField(lpResult, 0, COPTMEX_RESULT_DUAL, msol.dual);

    // 'qcslack'
    if (csol.nQConstr > 0) {
      bxAddField(lpResult, COPTMEX_RESULT_QCSLACK);
      bxSetField(lpResult, 0, COPTMEX_RESULT_QCSLACK, msol.qcslack);
    }

    if (csol.nPSD > 0) {
      // 'psdx'
      bxAddField(lpResult, COPTMEX_RESULT_PSDX);
      bxSetField(lpResult, 0, COPTMEX_RESULT_PSDX, msol.psdcolvalue);

      // 'psdrc'
      bxAddField(lpResult, COPTMEX_RESULT_PSDRC);
      bxSetField(lpResult, 0, COPTMEX_RESULT_PSDRC, msol.psdcoldual);
    }

    if (csol.nPSDConstr > 0) {
      // 'psdslack'
      bxAddField(lpResult, COPTMEX_RESULT_PSDSLACK);
      bxSetField(lpResult, 0, COPTMEX_RESULT_PSDSLACK, msol.psdrowslack);

      // 'psdpi'
      bxAddField(lpResult, COPTMEX_RESULT_PSDPI);
      bxSetField(lpResult, 0, COPTMEX_RESULT_PSDPI, msol.psdrowdual);
    }
  }

  if (msol.farkas != NULL) {
    // 'dualfarkas'
    bxAddField(lpResult, COPTMEX_RESULT_DUALFARKAS);
    bxSetField(lpResult, 0, COPTMEX_RESULT_DUALFARKAS, msol.farkas);
  }

  if (msol.ray != NULL) {
    // 'primalray'
    bxAddField(lpResult, COPTMEX_RESULT_PRIMALRAY);
    bxSetField(lpResult, 0, COPTMEX_RESULT_PRIMALRAY, msol.ray);
  }

  if (csol.hasBasis) {
    // 'varbasis'
    bxAddField(lpResult, COPTMEX_RESULT_VARBASIS);
    bxSetField(lpResult, 0, COPTMEX_RESULT_VARBASIS, msol.varbasis);
    // 'constrbasis'
    bxAddField(lpResult, COPTMEX_RESULT_CONBASIS);
    bxSetField(lpResult, 0, COPTMEX_RESULT_CONBASIS, msol.constrbasis);
  }

  *out_lpresult = lpResult;

exit_cleanup:
  if (retcode != COPT_RETCODE_OK) {
    *out_lpresult = NULL;
  }

  return retcode;
}

/* Extract MIP solution */
static int COPTMEX_getMipResult(copt_prob *prob, bxArray **out_mipresult) {
  int retcode = COPT_RETCODE_OK;
  bxArray *mipResult = NULL;
  coptmex_cmipsol csol;
  coptmex_mmipsol msol;

  COPTMEX_initCMipSol(&csol);
  COPTMEX_initMMipSol(&msol);

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &csol.nRow));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &csol.nCol));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASMIPSOL, &csol.hasMipSol));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_POOLSOLS, &csol.nSolPool));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_MIPSTATUS, &csol.nStatus));

  msol.status      = bxCreateString(COPTMEX_statusInt2Str(csol.nStatus));
  msol.simplexiter = bxCreateDoubleMatrix(1, 1, bxREAL);
  msol.nodecnt     = bxCreateDoubleMatrix(1, 1, bxREAL);
  msol.solvingtime = bxCreateDoubleMatrix(1, 1, bxREAL);
  if (!msol.status || !msol.simplexiter || !msol.nodecnt || !msol.solvingtime) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  if (csol.hasMipSol) {
    msol.bestgap = bxCreateDoubleMatrix(1, 1, bxREAL);
    msol.objval  = bxCreateDoubleMatrix(1, 1, bxREAL);
    msol.bestbnd = bxCreateDoubleMatrix(1, 1, bxREAL);
    msol.value   = bxCreateDoubleMatrix(csol.nCol, 1, bxREAL);
    if (!msol.bestgap || !msol.objval || !msol.bestbnd || !msol.value) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    csol.colValue = bxGetDoubles(msol.value);
  }

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_SIMPLEXITER, &csol.nSimplexIter));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_NODECNT, &csol.nNodeCnt));
  COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_SOLVINGTIME, &csol.dSolvingTime));

  if (csol.hasMipSol) {
    COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_BESTGAP, &csol.dBestGap));
    COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_BESTOBJ, &csol.dObjVal));
    COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_BESTBND, &csol.dBestBnd));
    COPTMEX_CALL(COPT_GetSolution(prob, csol.colValue));
  }

  if (csol.nSolPool > 0) {
    const char *solpoolfields[] = {COPTMEX_RESULT_POOLOBJ,
                                   COPTMEX_RESULT_POOLXN};
    msol.solpool = bxCreateStructMatrix(csol.nSolPool, 1, 2, solpoolfields);
    if (!msol.solpool) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    for (int i = 0; i < csol.nSolPool; ++i) {
      bxArray *poolobjval = bxCreateDoubleMatrix(1, 1, bxREAL);
      bxArray *poolxn = bxCreateDoubleMatrix(csol.nCol, 1, bxREAL);
      if (!poolobjval || !poolxn) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      double dPoolObjVal = +COPT_INFINITY;
      double *dPoolColValue = bxGetDoubles(poolxn);

      COPTMEX_CALL(COPT_GetPoolObjVal(prob, i, &dPoolObjVal));
      COPTMEX_CALL(COPT_GetPoolSolution(prob, i, csol.nCol, NULL, dPoolColValue));

      *bxGetDoubles(poolobjval) = dPoolObjVal;

      bxSetField(msol.solpool, i, COPTMEX_RESULT_POOLOBJ, poolobjval);
      bxSetField(msol.solpool, i, COPTMEX_RESULT_POOLXN, poolxn);
    }
  }

  *bxGetDoubles(msol.simplexiter) = csol.nSimplexIter;
  *bxGetDoubles(msol.nodecnt)     = csol.nNodeCnt;
  *bxGetDoubles(msol.solvingtime) = csol.dSolvingTime;

  if (csol.hasMipSol) {
    *bxGetDoubles(msol.bestgap) = csol.dBestGap;
    *bxGetDoubles(msol.objval)  = csol.dObjVal; 
    *bxGetDoubles(msol.bestbnd) = csol.dBestBnd;
  }

  mipResult = bxCreateStructMatrix(1, 1, 0, NULL);
  if (!mipResult) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // 'status'
  bxAddField(mipResult, COPTMEX_RESULT_STATUS);
  bxSetField(mipResult, 0, COPTMEX_RESULT_STATUS, msol.status);
  // 'simplexiter'
  bxAddField(mipResult, COPTMEX_RESULT_SIMITER);
  bxSetField(mipResult, 0, COPTMEX_RESULT_SIMITER, msol.simplexiter);
  // 'nodecnt'
  bxAddField(mipResult, COPTMEX_RESULT_NODECNT);
  bxSetField(mipResult, 0, COPTMEX_RESULT_NODECNT, msol.nodecnt);
  // 'solvingtime'
  bxAddField(mipResult, COPTMEX_RESULT_SOLVETIME);
  bxSetField(mipResult, 0, COPTMEX_RESULT_SOLVETIME, msol.solvingtime);

  if (csol.hasMipSol) {
    // 'bestgap'
    bxAddField(mipResult, COPTMEX_RESULT_BESTGAP);
    bxSetField(mipResult, 0, COPTMEX_RESULT_BESTGAP, msol.bestgap);
    // 'objval'
    bxAddField(mipResult, COPTMEX_RESULT_OBJVAL);
    bxSetField(mipResult, 0, COPTMEX_RESULT_OBJVAL, msol.objval);
    // 'bestbnd'
    bxAddField(mipResult, COPTMEX_RESULT_BESTBND);
    bxSetField(mipResult, 0, COPTMEX_RESULT_BESTBND, msol.bestbnd);
    // 'x'
    bxAddField(mipResult, COPTMEX_RESULT_VALUE);
    bxSetField(mipResult, 0, COPTMEX_RESULT_VALUE, msol.value);
  }

  if (csol.nSolPool > 0) {
    // 'pool'
    bxAddField(mipResult, COPTMEX_RESULT_POOL);
    bxSetField(mipResult, 0, COPTMEX_RESULT_POOL, msol.solpool);
  }

  *out_mipresult = mipResult;

exit_cleanup:
  if (retcode != COPT_RETCODE_OK) {
    *out_mipresult = NULL;
  }

  return retcode;
}

/* Extract and save result */
int COPTMEX_getResult(copt_prob *prob, bxArray **out_result) {
  int retcode = 0;
  int isMip = 0;

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ISMIP, &isMip));
  if (isMip) {
    COPTMEX_CALL(COPTMEX_getMipResult(prob, out_result));
  } else {
    COPTMEX_CALL(COPTMEX_getLpResult(prob, out_result));
  }

exit_cleanup:
  return retcode;
}

/* Extract model data */
int COPTMEX_getModel(copt_prob *prob, int nfiles, const bxArray **in_files,
                     bxArray **out_model) {
  int retcode = COPT_RETCODE_OK;
  int hasInfoFile = 0;
  bxArray *retmodel = NULL;
  coptmex_cprob cprob;
  coptmex_mprob mprob;

  // Read model from file
  if (nfiles == 1 || nfiles == 2) {
    COPTMEX_CALL(COPTMEX_readModel(prob, in_files[0]));
    if (nfiles == 2) {
      COPTMEX_CALL(COPTMEX_readInfo(prob, in_files[1]));
      hasInfoFile = 1;
    }
  }

  COPTMEX_initCProb(&cprob);
  COPTMEX_initMProb(&mprob);

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &cprob.nRow));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &cprob.nCol));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ELEMS, &cprob.nElem));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_QELEMS, &cprob.nQElem));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_SOSS, &cprob.nSos));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_INDICATORS, &cprob.nIndicator));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_CONES, &cprob.nCone));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_QCONSTRS, &cprob.nQConstr));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_OBJSENSE, &cprob.nObjSen));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASBASIS, &cprob.hasBasis));

  mprob.objsen = bxCreateString(COPTMEX_objsenInt2Str(cprob.nObjSen));
  mprob.objcon = bxCreateDoubleMatrix(1, 1, bxREAL);
  if (!mprob.objsen || !mprob.objcon) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_OBJCONST, &cprob.dObjConst));
  *bxGetDoubles(mprob.objcon) = cprob.dObjConst;

  if (cprob.nRow > 0 && cprob.nCol > 0) {
    int nRealElem = 0;
    if (cprob.nElem == 0) {
      nRealElem = 1;
    } else {
      nRealElem = cprob.nElem;
    }
    
    mprob.A = bxCreateSparse(cprob.nRow, cprob.nCol, nRealElem, bxREAL);
    if (!mprob.A) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    cprob.colMatBeg = (int *) bxCalloc(cprob.nCol + 1, sizeof(int));
    cprob.colMatIdx = (int *) bxCalloc(nRealElem, sizeof(int));
    if (!cprob.colMatBeg || !cprob.colMatIdx) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
    cprob.colMatElem = bxGetDoubles(mprob.A);

    if (cprob.nElem > 0) {
      COPTMEX_CALL(COPT_GetCols(prob, cprob.nCol, NULL, cprob.colMatBeg, NULL,
                   cprob.colMatIdx, cprob.colMatElem, cprob.nElem, NULL));
    }
    
    mwIndex *colMatBeg_data = bxGetJc(mprob.A);
    mwIndex *colMatIdx_data = bxGetIr(mprob.A);

    for (int i = 0; i < cprob.nCol + 1; ++i) {
      colMatBeg_data[i] = cprob.colMatBeg[i];
    }
    for (int i = 0; i < nRealElem; ++i) {
      colMatIdx_data[i] = cprob.colMatIdx[i];
    }

    bxFree(cprob.colMatBeg);
    bxFree(cprob.colMatIdx);
  }

  if (cprob.nCol > 0) {
    mprob.obj = bxCreateDoubleMatrix(cprob.nCol, 1, bxREAL);
    mprob.lb  = bxCreateDoubleMatrix(cprob.nCol, 1, bxREAL);
    mprob.ub  = bxCreateDoubleMatrix(cprob.nCol, 1, bxREAL);
    mprob.varnames = bxCreateCellMatrix(cprob.nCol, 1);
    if (!mprob.obj || !mprob.lb || !mprob.ub || !mprob.varnames) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    cprob.colCost  = bxGetDoubles(mprob.obj);
    cprob.colLower = bxGetDoubles(mprob.lb);
    cprob.colUpper = bxGetDoubles(mprob.ub);

    COPTMEX_CALL(COPT_GetColInfo(prob, COPT_DBLINFO_OBJ, cprob.nCol, NULL, cprob.colCost));
    COPTMEX_CALL(COPT_GetColInfo(prob, COPT_DBLINFO_LB, cprob.nCol, NULL, cprob.colLower));
    COPTMEX_CALL(COPT_GetColInfo(prob, COPT_DBLINFO_UB, cprob.nCol, NULL, cprob.colUpper));

    char *colType_c  = (char *) bxCalloc(2 * cprob.nCol, sizeof(char));
    char **colType_s = (char **) bxCalloc(cprob.nCol, sizeof(char *));
    if (!colType_c || !colType_s) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetColType(prob, cprob.nCol, NULL, colType_c));

    for (int i = cprob.nCol - 1; i >= 0; --i) {
      colType_c[2 * i] = colType_c[i];
      colType_c[2 * i + 1] = 0;
    }

    for (int i = 0; i < cprob.nCol; ++i) {
      colType_s[i] = colType_c + 2 * i;
    }

    mprob.vtype = bxCreateCharMatrixFromStrings(cprob.nCol, (const char **) colType_s);

    int colNameLen = 0;
    char **colNames = (char **) bxCalloc(cprob.nCol, sizeof(char *));
    if (!colNames) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    for (int i = 0; i < cprob.nCol; ++i) {
      COPTMEX_CALL(COPT_GetColName(prob, i, NULL, 0, &colNameLen));
      colNames[i] = (char *) bxCalloc(colNameLen, sizeof(char));
      if (!colNames[i]) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }
      COPTMEX_CALL(COPT_GetColName(prob, i, colNames[i], colNameLen, NULL));

      bxSetCell(mprob.varnames, i, bxCreateString(colNames[i]));
    }
  }

  if (cprob.nRow > 0) {
    mprob.lhs = bxCreateDoubleMatrix(cprob.nRow, 1, bxREAL);
    mprob.rhs = bxCreateDoubleMatrix(cprob.nRow, 1, bxREAL);
    mprob.constrnames = bxCreateCellMatrix(cprob.nRow, 1);
    if (!mprob.lhs || !mprob.rhs || !mprob.constrnames) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    cprob.rowLower = bxGetDoubles(mprob.lhs);
    cprob.rowUpper = bxGetDoubles(mprob.rhs);

    COPTMEX_CALL(COPT_GetRowInfo(prob, COPT_DBLINFO_LB, cprob.nRow, NULL, cprob.rowLower));
    COPTMEX_CALL(COPT_GetRowInfo(prob, COPT_DBLINFO_UB, cprob.nRow, NULL, cprob.rowUpper));
    
    int rowNameLen = 0;
    //TODO: rowSense
    char **rowNames = (char **) bxCalloc(cprob.nRow, sizeof(char *));
    if (!rowNames) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    for (int i = 0; i < cprob.nRow; ++i) {
      COPTMEX_CALL(COPT_GetRowName(prob, i, NULL, 0, &rowNameLen));
      rowNames[i] = (char *) bxCalloc(rowNameLen, sizeof(char));
      if (!rowNames[i]) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }
      COPTMEX_CALL(COPT_GetRowName(prob, i, rowNames[i], rowNameLen, NULL));

      bxSetCell(mprob.constrnames, i, bxCreateString(rowNames[i]));
    }
  }

  if (cprob.nSos > 0) {
    const char *sosfields[] = {COPTMEX_MODEL_SOSTYPE,
                               COPTMEX_MODEL_SOSVARS,
                               COPTMEX_MODEL_SOSWEIGHT};
    mprob.sos = bxCreateStructMatrix(cprob.nSos, 1, 3, sosfields);
    if (!mprob.sos) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetSOSs(prob, cprob.nSos, NULL, NULL, NULL, NULL, NULL,
                 NULL, 0, &cprob.nSosSize));
    
    cprob.sosType   = (int *) bxCalloc(cprob.nSos, sizeof(int));
    cprob.sosMatBeg = (int *) bxCalloc(cprob.nSos, sizeof(int));
    cprob.sosMatCnt = (int *) bxCalloc(cprob.nSos, sizeof(int));
    cprob.sosMatIdx = (int *) bxCalloc(cprob.nSosSize, sizeof(int));
    cprob.sosMatWt  = (double *) bxCalloc(cprob.nSosSize, sizeof(double));
    if (!cprob.sosType || !cprob.sosMatBeg || !cprob.sosMatCnt ||
        !cprob.sosMatIdx || !cprob.sosMatWt) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetSOSs(prob, cprob.nSos, NULL, cprob.sosType,
                 cprob.sosMatBeg, cprob.sosMatCnt, cprob.sosMatIdx,
                 cprob.sosMatWt, cprob.nSosSize, NULL));
    
    for (int i = 0; i < cprob.nSos; ++i) {
      bxArray *sosType = bxCreateDoubleMatrix(1, 1, bxREAL);
      bxArray *sosIdx  = bxCreateDoubleMatrix(cprob.sosMatCnt[i], 1, bxREAL);
      bxArray *sosWts  = bxCreateDoubleMatrix(cprob.sosMatCnt[i], 1, bxREAL);
      if (!sosType || !sosIdx || !sosWts) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      *bxGetDoubles(sosType) = cprob.sosType[i];

      int iSosElem = cprob.sosMatBeg[i];
      int iSosLast = cprob.sosMatCnt[i] + iSosElem;
      double *sosIdx_data = bxGetDoubles(sosIdx);
      double *sosWts_data = bxGetDoubles(sosWts);
      for (int iElem = 0; iElem < cprob.sosMatCnt[i]; ++iElem) {
        sosIdx_data[iElem] = cprob.sosMatIdx[iSosElem] + 1;
        sosWts_data[iElem] = cprob.sosMatWt[iSosElem];
        iSosElem++;
      }

      bxSetField(mprob.sos, i, COPTMEX_MODEL_SOSTYPE, sosType);
      bxSetField(mprob.sos, i, COPTMEX_MODEL_SOSVARS, sosIdx);
      bxSetField(mprob.sos, i, COPTMEX_MODEL_SOSWEIGHT, sosWts);
    }
  }

  if (cprob.nIndicator > 0) {
    const char *indicfields[] = {COPTMEX_MODEL_INDICBINVAR,
                                 COPTMEX_MODEL_INDICBINVAL,
                                 COPTMEX_MODEL_INDICROW,
                                 COPTMEX_MODEL_INDICSENSE,
                                 COPTMEX_MODEL_INDICRHS};
    mprob.indicator = bxCreateStructMatrix(cprob.nIndicator, 1, 5, indicfields);
    if (!mprob.indicator) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    for (int i = 0; i < cprob.nIndicator; ++i) {
      int rowElemCnt = 0;
      COPTMEX_CALL(COPT_GetIndicator(prob, i, NULL, NULL, NULL, NULL, NULL, 
                   NULL, NULL, 0, &rowElemCnt));

      bxArray *binVar = bxCreateDoubleMatrix(1, 1, bxREAL);
      bxArray *binVal = bxCreateDoubleMatrix(1, 1, bxREAL);
      bxArray *indicA = bxCreateSparse(cprob.nCol, 1, rowElemCnt, bxREAL);
      bxArray *indicSense = NULL;
      bxArray *indicRhs = bxCreateDoubleMatrix(1, 1, bxREAL);

      int binColIdx = 0;
      int binColVal = 0;
      int nRowMatCnt = 0;
      int *rowMatIdx = (int *) bxCalloc(rowElemCnt, sizeof(int));
      double *rowMatElem = (double *) bxCalloc(rowElemCnt, sizeof(double));
      char cRowSense[2] = {0};
      double dRowBound = 0;

      COPTMEX_CALL(COPT_GetIndicator(prob, i, &binColIdx, &binColVal,
                   &nRowMatCnt, rowMatIdx, rowMatElem, &cRowSense[0],
                   &dRowBound, rowElemCnt, NULL));

      *bxGetDoubles(binVar) = binColIdx + 1;
      *bxGetDoubles(binVal) = binColVal;

      mwIndex *jc_data = bxGetJc(indicA);
      mwIndex *ir_data = bxGetIr(indicA);
      double *val_data = bxGetDoubles(indicA);

      jc_data[0] = 0;
      jc_data[1] = rowElemCnt;
      for (int i = 0; i < rowElemCnt; ++i) {
        ir_data[i] = rowMatIdx[i];
        val_data[i] = rowMatElem[i];
      }

      indicSense = bxCreateString(cRowSense);
      *bxGetDoubles(indicRhs) = dRowBound;

      bxFree(rowMatIdx);
      bxFree(rowMatElem);

      bxSetField(mprob.indicator, i, COPTMEX_MODEL_INDICBINVAR, binVar);
      bxSetField(mprob.indicator, i, COPTMEX_MODEL_INDICBINVAL, binVal);
      bxSetField(mprob.indicator, i, COPTMEX_MODEL_INDICROW, indicA);
      bxSetField(mprob.indicator, i, COPTMEX_MODEL_INDICSENSE, indicSense);
      bxSetField(mprob.indicator, i, COPTMEX_MODEL_INDICRHS, indicRhs);
    }
  }

  if (cprob.nCone > 0) {
    const char *conefields[] = {COPTMEX_MODEL_CONETYPE,
                                COPTMEX_MODEL_CONEVARS};
    mprob.cone = bxCreateStructMatrix(cprob.nCone, 1, 2, conefields);
    if (!mprob.cone) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetCones(prob, cprob.nCone, NULL, NULL, NULL, NULL, NULL, 
                 0, &cprob.nConeSize));

    cprob.coneType = (int *) bxCalloc(cprob.nCone, sizeof(int));
    cprob.coneBeg  = (int *) bxCalloc(cprob.nCone, sizeof(int));
    cprob.coneCnt  = (int *) bxCalloc(cprob.nCone, sizeof(int));
    cprob.coneIdx  = (int *) bxCalloc(cprob.nConeSize, sizeof(int));
    if (!cprob.coneType || !cprob.coneBeg || !cprob.coneCnt || !cprob.coneIdx) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetCones(prob, cprob.nCone, NULL, cprob.coneType, 
                 cprob.coneBeg, cprob.coneCnt, cprob.coneIdx, cprob.nConeSize, 
                 NULL));

    for (int i = 0; i < cprob.nCone; ++i) {
      bxArray *coneType = bxCreateDoubleMatrix(1, 1, bxREAL);
      bxArray *coneIdx = bxCreateDoubleMatrix(cprob.coneCnt[i], 1, bxREAL);
      if (!coneType || !coneIdx) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      *bxGetDoubles(coneType) = cprob.coneType[i];

      int iConeElem = cprob.coneBeg[i];
      int iConeLast = cprob.coneCnt[i] + iConeElem;
      double *coneIdx_data = bxGetDoubles(coneIdx);
      for (int iElem = 0; iElem < cprob.coneCnt[i]; ++iElem) {
        coneIdx_data[iElem] = cprob.coneIdx[iConeElem] + 1;
        iConeElem++;
      }

      bxSetField(mprob.cone, i, COPTMEX_MODEL_CONETYPE, coneType);
      bxSetField(mprob.cone, i, COPTMEX_MODEL_CONEVARS, coneIdx);
    }
  }

  if (cprob.nQElem > 0) {
    mprob.qobj = bxCreateSparse(cprob.nCol, cprob.nCol, cprob.nQElem, bxREAL);
    if (!mprob.qobj) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    int *qObjRow = (int *) bxCalloc(cprob.nQElem, sizeof(int));
    int *qObjCol = (int *) bxCalloc(cprob.nQElem, sizeof(int));
    double *qObjElem = (double *) bxCalloc(cprob.nQElem, sizeof(double));
    if (!qObjRow || !qObjCol || !qObjElem) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetQuadObj(prob, NULL, qObjRow, qObjCol, qObjElem));
    COPTMEX_CALL(COPTMEX_coo2csc(cprob.nQElem, qObjRow, qObjCol, qObjElem, mprob.qobj));

    bxFree(qObjRow);
    bxFree(qObjCol);
    bxFree(qObjElem);
  }

  if (cprob.nQConstr > 0) {
    const char *qconstrfields[] = {COPTMEX_MODEL_QCROW,
                                   COPTMEX_MODEL_QCCOL,
                                   COPTMEX_MODEL_QCVAL,
                                   COPTMEX_MODEL_QCLINEAR,
                                   COPTMEX_MODEL_QCSENSE,
                                   COPTMEX_MODEL_QCRHS,
                                   COPTMEX_MODEL_QCNAME};
    mprob.quadcon = bxCreateStructMatrix(cprob.nQConstr, 1, 7, qconstrfields);
    if (!mprob.quadcon) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    for (int i = 0; i < cprob.nQConstr; ++i) {
      int nQMatElem = 0;
      int nQRowElem = 0;
      COPTMEX_CALL(COPT_GetQConstr(prob, i, NULL, NULL, NULL, 0, &nQMatElem, NULL,
                   NULL, NULL, NULL, 0, &nQRowElem));

      bxArray *QcRow = bxCreateDoubleMatrix(nQMatElem, 1, bxREAL);
      bxArray *QcCol = bxCreateDoubleMatrix(nQMatElem, 1, bxREAL);
      bxArray *QcVal = bxCreateDoubleMatrix(nQMatElem, 1, bxREAL);
      bxArray *QcLinear = bxCreateSparse(cprob.nCol, 1, nQRowElem, bxREAL);
      bxArray *QcSense = NULL;
      bxArray *QcRhs = bxCreateDoubleMatrix(1, 1, bxREAL);
      bxArray *QcName = NULL;

      int *qMatRow = (int *) bxCalloc(nQMatElem, sizeof(int));
      int *qMatCol = (int *) bxCalloc(nQMatElem, sizeof(int));
      double *qMatElem = bxGetDoubles(QcVal);
      int *qRowMatIdx = (int *) bxCalloc(nQRowElem, sizeof(int));
      double *qRowMatElem = bxGetDoubles(QcLinear);
      char qRowSense[2];
      double qRowBound = 0.0;
      char *qRowName = NULL;

      COPTMEX_CALL(COPT_GetQConstr(prob, i, qMatRow, qMatCol, qMatElem, nQMatElem, NULL,
                   qRowMatIdx, qRowMatElem, qRowSense, &qRowBound, nQRowElem, NULL));

      double *qMatRow_data = bxGetDoubles(QcRow);
      double *qMatCol_data = bxGetDoubles(QcCol);
      for (int i = 0; i < nQMatElem; ++i) {
        qMatRow_data[i] = qMatRow[i];
        qMatCol_data[i] = qMatCol[i];
      }

      mwIndex *jc_data = bxGetJc(QcLinear);
      mwIndex *ir_data = bxGetIr(QcLinear);

      jc_data[0] = 0;
      jc_data[1] = nQRowElem;
      for (int i = 0; i < nQRowElem; ++i) {
        ir_data[i] = qRowMatIdx[i];
      }

      QcSense = bxCreateString(qRowSense);
      *bxGetDoubles(QcRhs) = qRowBound;

      int nQcNameSize = 0;
      COPTMEX_CALL(COPT_GetQConstrName(prob, i, NULL, 0, &nQcNameSize));

      qRowName = (char *) bxCalloc(nQcNameSize + 1, sizeof(char));
      if (!qRowName) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      COPTMEX_CALL(COPT_GetQConstrName(prob, i, qRowName, nQcNameSize, NULL));
      QcName = bxCreateString(qRowName);

      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCROW, QcRow);
      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCCOL, QcCol);
      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCVAL, QcVal);
      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCLINEAR, QcLinear);
      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCSENSE, QcSense);
      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCRHS, QcRhs);
      bxSetField(mprob.quadcon, i, COPTMEX_MODEL_QCNAME, QcName);

      bxFree(qMatRow);
      bxFree(qMatCol);
      bxFree(qRowMatIdx);
      bxFree(qRowName);
    }
  }

  if (hasInfoFile && cprob.hasBasis) {
    mprob.varbasis = bxCreateDoubleMatrix(cprob.nCol, 1, bxREAL);
    mprob.constrbasis = bxCreateDoubleMatrix(cprob.nRow, 1, bxREAL);
    if (!mprob.varbasis || !mprob.constrbasis) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    cprob.colBasis = (int *) bxCalloc(cprob.nCol, sizeof(int));
    cprob.rowBasis = (int *) bxCalloc(cprob.nRow, sizeof(int));
    if (!cprob.colBasis || !cprob.rowBasis) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_CALL(COPT_GetColBasis(prob, cprob.nCol, NULL, cprob.colBasis));
    COPTMEX_CALL(COPT_GetRowBasis(prob, cprob.nRow, NULL, cprob.rowBasis));

    double *colBasis_data = bxGetDoubles(mprob.varbasis);
    double *rowBasis_data = bxGetDoubles(mprob.constrbasis);
    for (int i = 0; i < cprob.nCol; ++i) {
      colBasis_data[i] = cprob.colBasis[i];
    }
    for (int i = 0; i < cprob.nRow; ++i) {
      rowBasis_data[i] = cprob.rowBasis[i];
    }

    bxFree(cprob.colBasis);
    bxFree(cprob.rowBasis);
  }

  retmodel = bxCreateStructMatrix(1, 1, 0, NULL);
  if (!retmodel) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // 'objsen'
  bxAddField(retmodel, COPTMEX_MODEL_OBJSEN);
  bxSetField(retmodel, 0, COPTMEX_MODEL_OBJSEN, mprob.objsen);
  // 'objcon'
  bxAddField(retmodel, COPTMEX_MODEL_OBJCON);
  bxSetField(retmodel, 0, COPTMEX_MODEL_OBJCON, mprob.objcon);

  if (mprob.A != NULL) {
    // 'A'
    bxAddField(retmodel, COPTMEX_MODEL_A);
    bxSetField(retmodel, 0, COPTMEX_MODEL_A, mprob.A);
  }

  if (cprob.nCol > 0) {
    // 'obj'
    bxAddField(retmodel, COPTMEX_MODEL_OBJ);
    bxSetField(retmodel, 0, COPTMEX_MODEL_OBJ, mprob.obj);
    // 'lb'
    bxAddField(retmodel, COPTMEX_MODEL_LB);
    bxSetField(retmodel, 0, COPTMEX_MODEL_LB, mprob.lb);
    // 'ub'
    bxAddField(retmodel, COPTMEX_MODEL_UB);
    bxSetField(retmodel, 0, COPTMEX_MODEL_UB, mprob.ub);
    // 'vtype'
    bxAddField(retmodel, COPTMEX_MODEL_VTYPE);
    bxSetField(retmodel, 0, COPTMEX_MODEL_VTYPE, mprob.vtype);
    // 'varnames'
    bxAddField(retmodel, COPTMEX_MODEL_VARNAME);
    bxSetField(retmodel, 0, COPTMEX_MODEL_VARNAME, mprob.varnames);
  }

  if (cprob.nRow > 0) {
    //TODO: 'sense'
    // 'lhs'
    bxAddField(retmodel, COPTMEX_MODEL_LHS);
    bxSetField(retmodel, 0, COPTMEX_MODEL_LHS, mprob.lhs);
    // 'rhs'
    bxAddField(retmodel, COPTMEX_MODEL_RHS);
    bxSetField(retmodel, 0, COPTMEX_MODEL_RHS, mprob.rhs);
    // 'constrnames'
    bxAddField(retmodel, COPTMEX_MODEL_CONNAME);
    bxSetField(retmodel, 0, COPTMEX_MODEL_CONNAME, mprob.constrnames);
  }

  // 'sos'
  if (cprob.nSos > 0) {
    bxAddField(retmodel, COPTMEX_MODEL_SOS);
    bxSetField(retmodel, 0, COPTMEX_MODEL_SOS, mprob.sos);
  }

  // 'indicator'
  if (cprob.nIndicator > 0) {
    bxAddField(retmodel, COPTMEX_MODEL_INDICATOR);
    bxSetField(retmodel, 0, COPTMEX_MODEL_INDICATOR, mprob.indicator);
  }

  // 'cone'
  if (cprob.nCone > 0) {
    bxAddField(retmodel, COPTMEX_MODEL_CONE);
    bxSetField(retmodel, 0, COPTMEX_MODEL_CONE, mprob.cone);
  }

  // 'Q'
  if (cprob.nQElem > 0) {
    bxAddField(retmodel, COPTMEX_MODEL_QUADOBJ);
    bxSetField(retmodel, 0, COPTMEX_MODEL_QUADOBJ, mprob.qobj);
  }

  // 'quadcon'
  if (cprob.nQConstr > 0) {
    bxAddField(retmodel, COPTMEX_MODEL_QUADCON);
    bxSetField(retmodel, 0, COPTMEX_MODEL_QUADCON, mprob.quadcon);
  }

  if (hasInfoFile && cprob.hasBasis) {
    // 'varbasis'
    bxAddField(retmodel, COPTMEX_RESULT_VARBASIS);
    bxSetField(retmodel, 0, COPTMEX_RESULT_VARBASIS, mprob.varbasis);
    // 'constrbasis'
    bxAddField(retmodel, COPTMEX_RESULT_CONBASIS);
    bxSetField(retmodel, 0, COPTMEX_RESULT_CONBASIS, mprob.constrbasis);
  }

  *out_model = retmodel;

exit_cleanup:
  if (retcode != COPT_RETCODE_OK) {
    *out_model = NULL;
  }

  return retcode;
}

/* Load parameters to problem */
int COPTMEX_setParam(copt_prob *prob, const bxArray *in_param) {
  int retcode = 0;
  char msgbuf[COPT_BUFFSIZE];
  
  int islogging = 1;
  bxArray *logging = NULL;
  for (int i = bxGetNumberOfFields(in_param) - 1; i >= 0; --i) {
    const char *loggingname = bxGetFieldNameByNumber(in_param, i);
    if (mystrcmp(loggingname, COPT_INTPARAM_LOGGING) == 0) {
      logging = bxGetField(in_param, 0, loggingname);
      if (!bxIsScalar(logging) || bxIsChar(logging)) {
        snprintf(msgbuf, COPT_BUFFSIZE, "parameter.%s", COPT_INTPARAM_LOGGING);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }

      islogging = (int) bxGetScalar(logging);
      break;
    }
  }

  if (islogging == 1) {
    COPTMEX_CALL(COPTMEX_dispBanner());
  }
  
  if (logging != NULL) {
    COPTMEX_CALL(COPT_SetIntParam(prob, COPT_INTPARAM_LOGGING, islogging));
  }

  for (int i = 0; i < bxGetNumberOfFields(in_param); ++i) {
    int partype = -1;
    const char *parname = bxGetFieldNameByNumber(in_param, i);
    bxArray *pararray = bxGetField(in_param, 0, parname);

    if (mystrcmp(parname, COPT_INTPARAM_LOGGING) == 0) {
      continue;
    }

    COPTMEX_CALL(COPT_SearchParamAttr(prob, parname, &partype));
    if (partype != 0 && partype != 1) {
      snprintf(msgbuf, COPT_BUFFSIZE, "parameter.%s", parname);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NAME, msgbuf);
      goto exit_cleanup;
    }
    
    if (!bxIsScalar(pararray) || bxIsChar(pararray)) {
      snprintf(msgbuf, COPT_BUFFSIZE, "parameter.%s", parname);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }

    if (partype == 0) {
      COPTMEX_CALL(COPT_SetDblParam(prob, parname, bxGetScalar(pararray)));
    } else if (partype == 1) {
      COPTMEX_CALL(COPT_SetIntParam(prob, parname, (int) bxGetScalar(pararray)));
    }
  }

exit_cleanup:
  return retcode;
}

static char *COPTMEX_getFileExt(const char *filename) {
  char *tmpfilename = NULL;
  char *lastdot = NULL;
  int lenfile = strlen(filename);

  tmpfilename = (char *) bxCalloc(lenfile + 1, sizeof(char));
  memcpy(tmpfilename, filename, lenfile + 1);

  lastdot = strrchr(tmpfilename, '.');
  if (!lastdot || lastdot == tmpfilename) {
    return "";
  }

  if (lastdot != NULL) {
    if (strcmp(lastdot, ".gz") == 0) {
      *lastdot = '\0';

      lastdot = strrchr(tmpfilename, '.');
      if (!lastdot || lastdot == tmpfilename) {
        return "";
      }
    }
  }

  return lastdot + 1;
}

/* Read optional information from file */
int COPTMEX_readInfo(copt_prob *prob, const bxArray *in_info) {
  int retcode = 0;
  char *filename = NULL;
  char *fileext = NULL;

  COPTMEX_CALL(COPTMEX_getString(in_info, &filename));

  fileext = COPTMEX_getFileExt(filename);
  if (strcmp(fileext, "bas") == 0) {
    COPTMEX_CALL(COPT_ReadBasis(prob, filename));
  } else {
    retcode = COPT_RETCODE_INVALID;
  }

exit_cleanup:
  COPTMEX_freeString(&filename);
  return retcode;
}

/* Read model from file */
int COPTMEX_readModel(copt_prob *prob, const bxArray *in_model) {
  int retcode = 0;
  char *filename = NULL;
  char *fileext = NULL;

  COPTMEX_CALL(COPTMEX_getString(in_model, &filename));

  fileext = COPTMEX_getFileExt(filename);
  if (strcmp(fileext, "mps") == 0) {
    COPTMEX_CALL(COPT_ReadMps(prob, filename));
  } else if (strcmp(fileext, "lp") == 0) {
    COPTMEX_CALL(COPT_ReadLp(prob, filename));
  } else if (strcmp(fileext, "bin") == 0) {
    COPTMEX_CALL(COPT_ReadBin(prob, filename));
  } else if (strcmp(fileext, "dat-s") == 0) {
    COPTMEX_CALL(COPT_ReadSDPA(prob, filename));
  } else if (strcmp(fileext, "cbf") == 0) {
    COPTMEX_CALL(COPT_ReadCbf(prob, filename));
  } else {
    retcode = COPT_RETCODE_INVALID;
  }

exit_cleanup:
  COPTMEX_freeString(&filename);
  return retcode;
}

/* Write model to file */
int COPTMEX_writeModel(copt_prob *prob, const bxArray *out_file) {
  int retcode = 0;
  char *filename = NULL;
  char *fileext = NULL;

  COPTMEX_CALL(COPTMEX_getString(out_file, &filename));

  fileext = COPTMEX_getFileExt(filename);
  if (strcmp(fileext, "mps") == 0) {
    COPTMEX_CALL(COPT_WriteMps(prob, filename));
  } else if (strcmp(fileext, "lp") == 0) {
    COPTMEX_CALL(COPT_WriteLp(prob, filename));
  } else if (strcmp(fileext, "bin") == 0) {
    COPTMEX_CALL(COPT_WriteBin(prob, filename));
  } else if (strcmp(fileext, "cbf") == 0) {
    COPTMEX_CALL(COPT_WriteCbf(prob, filename));
  } else {
    retcode = COPT_RETCODE_INVALID;
  }

exit_cleanup:
  COPTMEX_freeString(&filename);
  return retcode;
}

/* Check if solve problem via cone data */
int COPTMEX_isConeModel(const bxArray *in_model) {
  int ifConeData = 0;
  bxArray *conedata = NULL;
  
  conedata = bxGetField(in_model, 0, COPTMEX_MODEL_CONEDATA);
  if (conedata != NULL) {
    if (COPTMEX_checkConeModel(conedata)) {
      ifConeData = 1;
    }
  }

  return ifConeData;
}

/* Solve cone problem with cone data */
int COPTMEX_solveConeModel(copt_prob *prob, const bxArray *in_model, bxArray **out_result, int ifRetResult) {
  int retcode = 0;
  coptmex_cconeprob cconeprob;
  coptmex_mconeprob mconeprob;
  bxArray *conedata = NULL;
  int *outRowMap = NULL;

  COPTMEX_initCConeProb(&cconeprob);
  COPTMEX_initMConeProb(&mconeprob);

  conedata = bxGetField(in_model, 0, COPTMEX_MODEL_CONEDATA);

  mconeprob.c = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_C);
  mconeprob.A = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_A);
  mconeprob.b = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_B);

  mconeprob.K = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_K);
  mconeprob.f = bxGetField(mconeprob.K, 0, COPTMEX_MODEL_CONEK_F);
  mconeprob.l = bxGetField(mconeprob.K, 0, COPTMEX_MODEL_CONEK_L);
  mconeprob.q = bxGetField(mconeprob.K, 0, COPTMEX_MODEL_CONEK_Q);
  mconeprob.r = bxGetField(mconeprob.K, 0, COPTMEX_MODEL_CONEK_R);
  mconeprob.s = bxGetField(mconeprob.K, 0, COPTMEX_MODEL_CONEK_S);

  mconeprob.objsen = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_OBJSEN);
  mconeprob.objcon = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_OBJCON);
  mconeprob.Q      = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_Q);

  // 'objsen'
  if (mconeprob.objsen != NULL) {
    COPTMEX_CALL(COPTMEX_getObjsen(mconeprob.objsen, &cconeprob.nObjSense));
  }
  // 'objcon'
  if (mconeprob.objcon != NULL) {
    cconeprob.dObjConst = bxGetScalar(mconeprob.objcon);
  }
  // 'Q'
  if (mconeprob.Q != NULL) {
    cconeprob.nQObjElem = bxGetNzmax(mconeprob.Q);

    cconeprob.qObjRow = (int *) bxCalloc(cconeprob.nQObjElem, sizeof(int));
    cconeprob.qObjCol = (int *) bxCalloc(cconeprob.nQObjElem, sizeof(int));
    cconeprob.qObjElem = (double *) bxCalloc(cconeprob.nQObjElem, sizeof(double));
    if (!cconeprob.qObjRow || !cconeprob.qObjCol || !cconeprob.qObjElem) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_csc2coo(mconeprob.Q, cconeprob.qObjRow, cconeprob.qObjCol, cconeprob.qObjElem);
  }

  // 'c'
  if (mconeprob.c != NULL) {
    cconeprob.colObj = bxGetDoubles(mconeprob.c);
  }
  // 'A'
  if (mconeprob.A != NULL) {
    cconeprob.nRow       = bxGetM(mconeprob.A);
    cconeprob.nCol       = bxGetN(mconeprob.A);
    cconeprob.nElem      = bxGetNzmax(mconeprob.A);
    cconeprob.colMatBeg  = (int *) bxCalloc(cconeprob.nCol + 1, sizeof(int));
    cconeprob.colMatIdx  = (int *) bxCalloc(cconeprob.nElem, sizeof(int));
    if (!cconeprob.colMatBeg || !cconeprob.colMatIdx) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    mwIndex *colMatBeg_data = bxGetJc(mconeprob.A);
    mwIndex *colMatIdx_data = bxGetIr(mconeprob.A);
    for (int i = 0; i < cconeprob.nCol + 1; ++i) {
      cconeprob.colMatBeg[i] = (int) colMatBeg_data[i];
    }
    for (int i = 0; i < cconeprob.nElem; ++i) {
      cconeprob.colMatIdx[i] = (int) colMatIdx_data[i];
    }

    cconeprob.colMatElem = bxGetDoubles(mconeprob.A);
  }
  // 'b'
  if (mconeprob.b != NULL) {
    cconeprob.rowRhs = bxGetDoubles(mconeprob.b);
  }

  // 'K'
  if (mconeprob.K != NULL) {
    // 'f'
    if (mconeprob.f != NULL) {
      cconeprob.nFree = (int) bxGetScalar(mconeprob.f);
    }
    // 'l'
    if (mconeprob.l != NULL) {
      cconeprob.nPositive = (int) bxGetScalar(mconeprob.l);
    }
    // 'q'
    if (mconeprob.q != NULL) {
      int nCone = bxGetNumberOfElements(mconeprob.q);
      double *coneDim_data = bxGetDoubles(mconeprob.q);

      if (nCone > 1 || (nCone == 1 && coneDim_data[0] > 0)) {
        cconeprob.nCone = nCone;
        cconeprob.coneDim = (int *) bxCalloc(cconeprob.nCone, sizeof(int));

        for (int i = 0; i < cconeprob.nCone; ++i) {
          cconeprob.coneDim[i] = (int) coneDim_data[i];
        }
      }
    }
    // 'r'
    if (mconeprob.r != NULL) {
      int nRotateCone = bxGetNumberOfElements(mconeprob.r);
      double *rotateConeDim_data = bxGetDoubles(mconeprob.r);

      if (nRotateCone > 1 || (nRotateCone == 1 && rotateConeDim_data[0] > 0)) {
        cconeprob.nRotateCone = nRotateCone;
        cconeprob.rotateConeDim = (int *) bxCalloc(cconeprob.nRotateCone, sizeof(int));

        for (int i = 0; i < cconeprob.nRotateCone; ++i) {
          cconeprob.rotateConeDim[i] = (int) rotateConeDim_data[i];
        }
      }
    }
    // 's'
    if (mconeprob.s != NULL) {
      int nPSD = bxGetNumberOfElements(mconeprob.s);
      double *psdDim_data = bxGetDoubles(mconeprob.s);

      if (nPSD > 1 || (nPSD == 1 && psdDim_data[0] > 0)) {
        cconeprob.nPSD = nPSD;
        cconeprob.psdDim = (int *) bxCalloc(cconeprob.nPSD, sizeof(int));

        for (int i = 0; i < cconeprob.nPSD; ++i) {
          cconeprob.psdDim[i] = (int) psdDim_data[i];
        }
      }
    }
  }

  outRowMap = (int *) bxCalloc(cconeprob.nRow, sizeof(int));
  if (!outRowMap) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // Load cone problem data
  COPTMEX_CALL(COPT_LoadConeProb(prob, cconeprob.nCol, cconeprob.nRow,
    cconeprob.nFree, cconeprob.nPositive, 0, cconeprob.nCone, cconeprob.nRotateCone,
    0, 0, 0, 0, cconeprob.nPSD, cconeprob.nQObjElem, cconeprob.nObjSense, cconeprob.dObjConst,
    cconeprob.colObj, cconeprob.qObjRow, cconeprob.qObjCol, cconeprob.qObjElem,
    cconeprob.colMatBeg, NULL, cconeprob.colMatIdx, cconeprob.colMatElem, cconeprob.rowRhs,
    NULL, NULL, cconeprob.coneDim, cconeprob.rotateConeDim, NULL, NULL, NULL, NULL,
    cconeprob.psdDim, NULL, NULL, NULL, NULL, outRowMap));

  COPTMEX_CALL(COPT_Solve(prob));

  if (ifRetResult == 1) {
    COPTMEX_CALL(COPTMEX_getResult(prob, out_result));

    if (*out_result != NULL) {
      bxArray *rowMap = bxCreateDoubleMatrix(cconeprob.nRow, 1, bxREAL);
      if (rowMap == NULL) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      double *rowMap_data = bxGetDoubles(rowMap);
      for (int i = 0; i < cconeprob.nRow; ++i) {
        rowMap_data[i] = outRowMap[i];
      }

      bxAddField(*out_result, "rowmap");
      bxSetField(*out_result, 0, "rowmap", rowMap);
    }
  }

exit_cleanup:
  if (outRowMap != NULL) {
    bxFree(outRowMap);
  }

  if (cconeprob.qObjRow != NULL) {
    bxFree(cconeprob.qObjRow);
  }
  if (cconeprob.qObjCol != NULL) {
    bxFree(cconeprob.qObjCol);
  }
  if (cconeprob.qObjElem != NULL) {
    bxFree(cconeprob.qObjElem);
  }

  if (cconeprob.coneDim != NULL) {
    bxFree(cconeprob.coneDim);
  }
  if (cconeprob.rotateConeDim != NULL) {
    bxFree(cconeprob.rotateConeDim);
  }
  if (cconeprob.psdDim != NULL) {
    bxFree(cconeprob.psdDim);
  }

  if (cconeprob.colMatBeg != NULL) {
    bxFree(cconeprob.colMatBeg);
  }
  if (cconeprob.colMatIdx != NULL) {
    bxFree(cconeprob.colMatIdx);
  }
  return retcode;
}

/* Extract and load data to problem */
int COPTMEX_loadModel(copt_prob *prob, const bxArray *in_model) {
  int retcode = 0;
  coptmex_cprob cprob;
  coptmex_mprob mprob;

  COPTMEX_initCProb(&cprob);
  COPTMEX_initMProb(&mprob);

  mprob.objsen      = bxGetField(in_model, 0, COPTMEX_MODEL_OBJSEN);
  mprob.objcon      = bxGetField(in_model, 0, COPTMEX_MODEL_OBJCON);
  mprob.A           = bxGetField(in_model, 0, COPTMEX_MODEL_A);
  mprob.obj         = bxGetField(in_model, 0, COPTMEX_MODEL_OBJ);
  mprob.lb          = bxGetField(in_model, 0, COPTMEX_MODEL_LB);
  mprob.ub          = bxGetField(in_model, 0, COPTMEX_MODEL_UB);
  mprob.vtype       = bxGetField(in_model, 0, COPTMEX_MODEL_VTYPE);
  mprob.varnames    = bxGetField(in_model, 0, COPTMEX_MODEL_VARNAME);
  mprob.sense       = bxGetField(in_model, 0, COPTMEX_MODEL_SENSE);
  mprob.lhs         = bxGetField(in_model, 0, COPTMEX_MODEL_LHS);
  mprob.rhs         = bxGetField(in_model, 0, COPTMEX_MODEL_RHS);
  mprob.constrnames = bxGetField(in_model, 0, COPTMEX_MODEL_CONNAME);

  mprob.sos         = bxGetField(in_model, 0, COPTMEX_MODEL_SOS);
  mprob.indicator   = bxGetField(in_model, 0, COPTMEX_MODEL_INDICATOR);
  mprob.cone        = bxGetField(in_model, 0, COPTMEX_MODEL_CONE);

  mprob.qobj        = bxGetField(in_model, 0, COPTMEX_MODEL_QUADOBJ);
  mprob.quadcon     = bxGetField(in_model, 0, COPTMEX_MODEL_QUADCON);

  mprob.varbasis    = bxGetField(in_model, 0, COPTMEX_RESULT_VARBASIS);
  mprob.constrbasis = bxGetField(in_model, 0, COPTMEX_RESULT_CONBASIS);

  mprob.value       = bxGetField(in_model, 0, COPTMEX_RESULT_VALUE);
  mprob.slack       = bxGetField(in_model, 0, COPTMEX_RESULT_SLACK);
  mprob.dual        = bxGetField(in_model, 0, COPTMEX_RESULT_DUAL);
  mprob.redcost     = bxGetField(in_model, 0, COPTMEX_RESULT_REDCOST);

  mprob.mipstart    = bxGetField(in_model, 0, COPTMEX_ADVINFO_MIPSTART);

  if (COPTMEX_checkModel(&mprob) == 0) {
    goto exit_cleanup;
  }

  // 'objsen'
  if (mprob.objsen != NULL) {
    COPTMEX_CALL(COPTMEX_getObjsen(mprob.objsen, &cprob.nObjSen));
  }
  // 'objcon'
  if (mprob.objcon != NULL) {
    cprob.dObjConst = bxGetScalar(mprob.objcon);
  }
  // 'A'
  if (mprob.A != NULL) {
    cprob.nRow       = bxGetM(mprob.A);
    cprob.nCol       = bxGetN(mprob.A);
    cprob.nElem      = bxGetNzmax(mprob.A);
    cprob.colMatBeg  = (int *) bxCalloc(cprob.nCol + 1, sizeof(int));
    cprob.colMatIdx  = (int *) bxCalloc(cprob.nElem, sizeof(int));
    if (!cprob.colMatBeg || !cprob.colMatIdx) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    mwIndex *colMatBeg_data = bxGetJc(mprob.A);
    mwIndex *colMatIdx_data = bxGetIr(mprob.A);
    for (int i = 0; i < cprob.nCol + 1; ++i) {
      cprob.colMatBeg[i] = (int) colMatBeg_data[i];
    }
    for (int i = 0; i < cprob.nElem; ++i) {
      cprob.colMatIdx[i] = (int) colMatIdx_data[i];
    }

    cprob.colMatElem = bxGetDoubles(mprob.A);
  }
  // 'obj'
  if (mprob.obj != NULL) {
    cprob.colCost = bxGetDoubles(mprob.obj);
  }
  // 'lb'
  if (mprob.lb != NULL) {
    cprob.colLower = bxGetDoubles(mprob.lb);
  }
  // 'ub'
  if (mprob.ub != NULL) {
    cprob.colUpper = bxGetDoubles(mprob.ub);
  }
  // 'vtype'
  if (mprob.vtype != NULL) {
    if (bxGetNumberOfElements(mprob.vtype) == cprob.nCol) {
      COPTMEX_CALL(COPTMEX_getString(mprob.vtype, &cprob.colType));
    } else {
      char *vtype = NULL;
      COPTMEX_CALL(COPTMEX_getString(mprob.vtype, &vtype));

      cprob.colType = (char *) bxCalloc(cprob.nCol + 1, sizeof(char));
      if (!cprob.colType) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }
      for (int i = 0; i < cprob.nCol; ++i) {
        cprob.colType[i] = vtype[0];
      }
    }
  }
  // 'varnames'
  if (mprob.varnames != NULL) {
    cprob.colNames = (char **) bxCalloc(cprob.nCol, sizeof(char *));
    for (int i = 0; i < cprob.nCol; ++i) {
      bxArray *nameCell = bxGetCell(mprob.varnames, i);
      COPTMEX_CALL(COPTMEX_getString(nameCell, &cprob.colNames[i]));
    }
  }
  // 'sense', 'lhs' and 'rhs'
  if (mprob.sense == NULL) {
    cprob.rowLower = bxGetDoubles(mprob.lhs);
    cprob.rowUpper = bxGetDoubles(mprob.rhs);
  } else {
    if (bxGetNumberOfElements(mprob.sense) == cprob.nRow) {
      COPTMEX_CALL(COPTMEX_getString(mprob.sense, &cprob.rowSense));
    } else {
      char *rsense = NULL;
      COPTMEX_CALL(COPTMEX_getString(mprob.sense, &rsense));

      cprob.rowSense = (char *) bxCalloc(cprob.nRow + 1, sizeof(char));
      if (!cprob.rowSense) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }
      for (int i = 0; i < cprob.nRow; ++i) {
        cprob.rowSense[i] = rsense[0];
      }
    }

    cprob.rowUpper = bxGetDoubles(mprob.rhs);
  }
  // 'constrnames'
  if (mprob.constrnames != NULL) {
    cprob.rowNames = (char **) bxCalloc(cprob.nRow, sizeof(char *));
    for (int i = 0; i < cprob.nRow; ++i) {
      bxArray *namecell = bxGetCell(mprob.constrnames, i);
      COPTMEX_CALL(COPTMEX_getString(namecell, &cprob.rowNames[i]));
    }
  }

  // Load problem data to COPT problem
  if (cprob.rowSense == NULL) {
    COPTMEX_CALL(COPT_LoadProb(prob, cprob.nCol, cprob.nRow,
                 cprob.nObjSen, cprob.dObjConst, cprob.colCost,
                 cprob.colMatBeg, NULL, cprob.colMatIdx, cprob.colMatElem,
                 cprob.colType, cprob.colLower, cprob.colUpper,
                 NULL, cprob.rowLower, cprob.rowUpper,
                 cprob.colNames, cprob.rowNames));
  } else {
    COPTMEX_CALL(COPT_LoadProb(prob, cprob.nCol, cprob.nRow,
                 cprob.nObjSen, cprob.dObjConst, cprob.colCost,
                 cprob.colMatBeg, NULL, cprob.colMatIdx, cprob.colMatElem,
                 cprob.colType, cprob.colLower, cprob.colUpper,
                 cprob.rowSense, cprob.rowUpper, NULL,
                 cprob.colNames, cprob.rowNames));
  }
  
  // Extract and load the optional SOS part
  if (mprob.sos != NULL) {
    for (int i = 0; i < bxGetNumberOfElements(mprob.sos); ++i) {
      bxArray *sostype_m = bxGetField(mprob.sos, i, COPTMEX_MODEL_SOSTYPE);
      bxArray *sosvars_m = bxGetField(mprob.sos, i, COPTMEX_MODEL_SOSVARS);
      bxArray *soswgts_m = bxGetField(mprob.sos, i, COPTMEX_MODEL_SOSWEIGHT);

      int sosType    = (int) bxGetScalar(sostype_m);
      int sosMatBeg  = 0;
      int sosMatCnt  = (int) bxGetNumberOfElements(sosvars_m);
      int *sosMatIdx = (int *) bxCalloc(sosMatCnt, sizeof(int));
      if (!sosMatIdx) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      double *sosvars_data = bxGetDoubles(sosvars_m);
      for (int i = 0; i < sosMatCnt; ++i) {
        sosMatIdx[i] = (int) sosvars_data[i] - 1;
      }

      double *sosMatWt = NULL;
      if (soswgts_m != NULL) {
        sosMatWt = bxGetDoubles(soswgts_m);
      }

      COPTMEX_CALL(COPT_AddSOSs(prob, 1, &sosType, &sosMatBeg, &sosMatCnt,
                   sosMatIdx, sosMatWt));
      
      bxFree(sosMatIdx);
    }
  }

  // Extract and load the optional indicator part
  if (mprob.indicator != NULL) {
    for (int i = 0; i < bxGetNumberOfElements(mprob.indicator); ++i) {
      bxArray *binVar = bxGetField(mprob.indicator, i, COPTMEX_MODEL_INDICBINVAR);
      bxArray *binVal = bxGetField(mprob.indicator, i, COPTMEX_MODEL_INDICBINVAL);
      bxArray *indicA = bxGetField(mprob.indicator, i, COPTMEX_MODEL_INDICROW);
      bxArray *rSense = bxGetField(mprob.indicator, i, COPTMEX_MODEL_INDICSENSE);
      bxArray *rowBnd = bxGetField(mprob.indicator, i, COPTMEX_MODEL_INDICRHS);

      int binColIdx  = (int) bxGetScalar(binVar) - 1;
      int binColVal  = (int) bxGetScalar(binVal);
      int nRowMatCnt = 0;
      int *rowMatIdx = NULL;
      double *rowMatElem = NULL;

      if (bxIsSparse(indicA)) {
        nRowMatCnt = bxGetNzmax(indicA);
        rowMatIdx  = (int *) bxCalloc(nRowMatCnt, sizeof(int));
        if (!rowMatIdx) {
          retcode = COPT_RETCODE_MEMORY;
          goto exit_cleanup;
        }
        mwIndex *rowMatIdx_data = bxGetIr(indicA);
        for (int i = 0; i < nRowMatCnt; ++i) {
          rowMatIdx[i] = rowMatIdx_data[i];
        }
        rowMatElem = bxGetDoubles(indicA);
      } else {
        double *rowMatElem_data = bxGetDoubles(indicA);
        for (int i = 0; i < bxGetNumberOfElements(indicA); ++i) {
          if (rowMatElem_data[i] != 0) {
            ++nRowMatCnt;
          }
        }
        rowMatIdx = (int *) bxCalloc(nRowMatCnt, sizeof(int));
        rowMatElem = (double *) bxCalloc(nRowMatCnt, sizeof(double));
        if (!rowMatIdx || !rowMatElem) {
          retcode = COPT_RETCODE_MEMORY;
          goto exit_cleanup;
        }
        for (int i = 0, iElem = 0; i < bxGetNumberOfElements(indicA); ++i) {
          if (rowMatElem_data[i] != 0) {
            rowMatIdx[iElem] = i;
            rowMatElem[iElem] = rowMatElem_data[i];
            iElem++;
          }
        }
      }

      char cRowSense[2];
      bxAsCStr(rSense, cRowSense, 2);
      double dRowBound = bxGetScalar(rowBnd);

      COPTMEX_CALL(COPT_AddIndicator(prob, binColIdx, binColVal, nRowMatCnt,
                   rowMatIdx, rowMatElem, cRowSense[0], dRowBound));

      bxFree(rowMatIdx);
      if (!bxIsSparse(indicA)) {
        bxFree(rowMatElem);
      }
    }
  }

  // Extract and load the optional cone part
  if (mprob.cone != NULL) {
    for (int i = 0; i < bxGetNumberOfElements(mprob.cone); ++i) {
      bxArray *conetype_m = bxGetField(mprob.cone, i, COPTMEX_MODEL_CONETYPE);
      bxArray *conevars_m = bxGetField(mprob.cone, i, COPTMEX_MODEL_CONEVARS);

      int coneType = (int) bxGetScalar(conetype_m);
      int coneBeg = 0;
      int coneCnt = (int) bxGetNumberOfElements(conevars_m);
      int *coneIdx = (int *) bxCalloc(coneCnt, sizeof(int));
      if (!coneIdx) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      double *conevars_data = bxGetDoubles(conevars_m);
      for (int i = 0; i < coneCnt; ++i) {
        coneIdx[i] = (int) conevars_data[i] - 1;
      }

      COPTMEX_CALL(COPT_AddCones(prob, 1, &coneType, &coneBeg, &coneCnt, coneIdx));

      bxFree(coneIdx);
    }
  }

  // Extract and load optional Q objective part
  if (mprob.qobj != NULL) {
    cprob.nQElem = bxGetNzmax(mprob.qobj);
    int *qObjRow = (int *) bxCalloc(cprob.nQElem, sizeof(int));
    int *qObjCol = (int *) bxCalloc(cprob.nQElem, sizeof(int));
    double *qObjElem = (double *) bxCalloc(cprob.nQElem, sizeof(double));
    if (!qObjRow || !qObjCol || !qObjElem) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }

    COPTMEX_csc2coo(mprob.qobj, qObjRow, qObjCol, qObjElem);
    COPTMEX_CALL(COPT_SetQuadObj(prob, cprob.nQElem, qObjRow, qObjCol, qObjElem));

    bxFree(qObjRow);
    bxFree(qObjCol);
    bxFree(qObjElem);
  }

  // Extract and load optional quadratic constraint part
  if (mprob.quadcon != NULL) {
    for (int i = 0; i < bxGetNumberOfElements(mprob.quadcon); ++i) {
      bxArray *QcMat = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCSPMAT);
      bxArray *QcRow = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCROW);
      bxArray *QcCol = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCCOL);
      bxArray *QcVal = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCVAL);
      bxArray *QcLinear = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCLINEAR);
      bxArray *QcSense = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCSENSE);
      bxArray *QcRhs = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCRHS);
      bxArray *QcName = bxGetField(mprob.quadcon, i, COPTMEX_MODEL_QCNAME);

      int nQMatElem = 0;
      int *qMatRow = NULL;
      int *qMatCol = NULL;
      double *qMatElem = NULL;
      int nQRowElem = 0;
      int *qRowMatIdx = NULL;
      double *qRowMatElem = NULL;
      char qRowSense[2];
      double qRowBound = 0.0;
      char qRowName[COPT_BUFFSIZE] = {0};

      if (QcMat != NULL) {
        nQMatElem = bxGetNzmax(QcMat);
        qMatRow = (int *) bxCalloc(nQMatElem, sizeof(int));
        qMatCol = (int *) bxCalloc(nQMatElem, sizeof(int));
        qMatElem = (double *) bxCalloc(nQMatElem, sizeof(double));
        if (!qMatRow || !qMatCol || !qMatElem) {
          retcode = COPT_RETCODE_MEMORY;
          goto exit_cleanup;
        }

        COPTMEX_csc2coo(QcMat, qMatRow, qMatCol, qMatElem);
      } else {
        if (QcRow != NULL && QcCol != NULL && QcVal != NULL) {
          nQMatElem = bxGetNumberOfElements(QcRow);
          qMatRow = (int *) bxCalloc(nQMatElem, sizeof(int));
          qMatCol = (int *) bxCalloc(nQMatElem, sizeof(int));
          if (!qMatRow || !qMatCol) {
            retcode = COPT_RETCODE_MEMORY;
            goto exit_cleanup;
          }

          double *qMatRow_data = bxGetDoubles(QcRow);
          double *qMatCol_data = bxGetDoubles(QcCol);
          qMatElem = bxGetDoubles(QcVal);
          for (int i = 0; i < nQMatElem; ++i) {
            qMatRow[i] = (int) qMatRow_data[i] - 1;
            qMatCol[i] = (int) qMatCol_data[i] - 1;
          }
        }
      }

      if (QcLinear != NULL) {
        if (bxIsSparse(QcLinear)) {
          double *qRowMatElem_data = bxGetDoubles(QcLinear);
          for (int i = 0; i < bxGetNzmax(QcLinear); ++i) {
            if (qRowMatElem_data[i] != 0.0) {
              ++nQRowElem;
            }
          }
          if (nQRowElem > 0) {
            qRowMatIdx = (int *) bxCalloc(nQRowElem, sizeof(int));
            qRowMatElem = (double *) bxCalloc(nQRowElem, sizeof(double));
            if (!qRowMatIdx || !qRowMatElem) {
              retcode = COPT_RETCODE_MEMORY;
              goto exit_cleanup;
            }
            mwIndex *qRowMatIdx_data = bxGetIr(QcLinear);
            for (int i = 0, iElem = 0; i < bxGetNzmax(QcLinear); ++i) {
              if (qRowMatElem_data[i] != 0.0) {
                qRowMatIdx[iElem] = (int) qRowMatIdx_data[i];
                qRowMatElem[iElem] = qRowMatElem_data[i];
                iElem++;
              }
            }
          }
        } else {
          double *qRowMatElem_data = bxGetDoubles(QcLinear);
          for (int i = 0; i < bxGetNumberOfElements(QcLinear); ++i) {
            if (qRowMatElem_data[i] != 0.0) {
              ++nQRowElem;
            }
          }
          if (nQRowElem > 0) {
            qRowMatIdx = (int *) bxCalloc(nQRowElem, sizeof(int));
            qRowMatElem = (double *) bxCalloc(nQRowElem, sizeof(double));
            if (!qRowMatIdx || !qRowMatElem) {
              retcode = COPT_RETCODE_MEMORY;
              goto exit_cleanup;
            }
            for (int i = 0, iElem = 0; i < bxGetNumberOfElements(QcLinear); ++i) {
              if (qRowMatElem_data[i] != 0.0) {
                qRowMatIdx[iElem] = i;
                qRowMatElem[iElem] = qRowMatElem_data[i];
                iElem++;
              }
            }
          }
        }
      }

      if (QcSense != NULL) {
        bxAsCStr(QcSense, qRowSense, 2);
      } else {
        qRowSense[0] = COPT_LESS_EQUAL;
      }
      qRowBound = bxGetScalar(QcRhs);
      if (QcName != NULL) {
        bxAsCStr(QcName, qRowName, COPT_BUFFSIZE);
      }

      COPTMEX_CALL(COPT_AddQConstr(prob, nQRowElem, qRowMatIdx, qRowMatElem, 
                   nQMatElem, qMatRow, qMatCol, qMatElem, qRowSense[0], qRowBound,
                   qRowName));

      bxFree(qMatRow);
      bxFree(qMatCol);
      if (QcMat != NULL) {
        bxFree(qMatElem);
      }

      if (nQRowElem > 0) {
        bxFree(qRowMatIdx);
        bxFree(qRowMatElem);
      }
    }
  }

  // Extract and load the optional advanced information
  if (mprob.varbasis != NULL && mprob.constrbasis != NULL) {
    cprob.colBasis = (int *) bxCalloc(cprob.nCol, sizeof(int));
    cprob.rowBasis = (int *) bxCalloc(cprob.nRow, sizeof(int));
    if (!cprob.colBasis || !cprob.rowBasis) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
    
    double *colBasis_data = bxGetDoubles(mprob.varbasis);
    double *rowBasis_data = bxGetDoubles(mprob.constrbasis);
    for (int i = 0; i < cprob.nCol; ++i) {
      cprob.colBasis[i] = (int) colBasis_data[i];
    }
    for (int i = 0; i < cprob.nRow; ++i) {
      cprob.rowBasis[i] = (int) rowBasis_data[i];
    }

    COPTMEX_CALL(COPT_SetBasis(prob, cprob.colBasis, cprob.rowBasis));

    bxFree(cprob.colBasis);
    bxFree(cprob.rowBasis);
  }

  if (mprob.value != NULL && mprob.slack != NULL && mprob.dual != NULL && mprob.redcost != NULL) {
    double *colValue = bxGetDoubles(mprob.value);
    double *colDual  = bxGetDoubles(mprob.redcost);
    double *rowSlack = bxGetDoubles(mprob.slack);
    double *rowDual  = bxGetDoubles(mprob.dual);

    COPTMEX_CALL(COPT_SetLpSolution(prob, colValue, rowSlack, rowDual, colDual));
  }

  if (mprob.mipstart != NULL) {
    int nRowCnt = 0;
    int *rowIdx = NULL;
    double *rowElem = NULL;

    if (bxIsSparse(mprob.mipstart)) {
      nRowCnt = bxGetNzmax(mprob.mipstart);
      rowIdx = (int *) bxCalloc(nRowCnt, sizeof(int));
      if (!rowIdx) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }
      mwIndex *rowIdx_data = bxGetIr(mprob.mipstart);
      for (int i = 0; i < nRowCnt; ++i) {
        rowIdx[i] = rowIdx_data[i];
      }
      rowElem = bxGetDoubles(mprob.mipstart);

      COPTMEX_CALL(COPT_AddMipStart(prob, nRowCnt, rowIdx, rowElem));

      bxFree(rowIdx);
    } else {
      double *rowElem_data = bxGetDoubles(mprob.mipstart);

      nRowCnt = bxGetNumberOfElements(mprob.mipstart);
      rowElem = (double *) bxCalloc(nRowCnt, sizeof(double));
      if (!rowElem) {
        retcode = COPT_RETCODE_MEMORY;
        goto exit_cleanup;
      }

      for (int i = 0; i < nRowCnt; ++i) {
        if (bxIsNaN(rowElem_data[i])) {
          rowElem[i] = COPT_UNDEFINED;
        } else {
          rowElem[i] = rowElem_data[i];
        }
      }

      COPTMEX_CALL(COPT_AddMipStart(prob, nRowCnt, NULL, rowElem));

      bxFree(rowElem); 
    }
  }

exit_cleanup:
  return retcode;
}

/* Extract IIS information */
static int COPTMEX_getIIS(copt_prob *prob, bxArray **out_iis) {
  int retcode = COPT_RETCODE_OK;
  int nRow = 0, nCol = 0, nSos = 0, nIndicator = 0;
  int hasIIS = 0;
  int isMinIIS = 0;
  bxArray *iisInfo = NULL;
  coptmex_ciisinfo ciisinfo;
  coptmex_miisinfo miisinfo;

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASIIS, &hasIIS));
  if (hasIIS == 0) {
    *out_iis = NULL;
    goto exit_cleanup;
  }

  COPTMEX_initCIISInfo(&ciisinfo);
  COPTMEX_initMIISInfo(&miisinfo);

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &nRow));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &nCol));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_SOSS, &nSos));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_INDICATORS, &nIndicator));

  miisinfo.isminiis = bxCreateDoubleMatrix(1, 1, bxREAL);
  if (!miisinfo.isminiis) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  if (nCol > 0) {
    miisinfo.varlb = bxCreateLogicalMatrix(nCol, 1);
    miisinfo.varub = bxCreateLogicalMatrix(nCol, 1);
    if (!miisinfo.varlb || !miisinfo.varub) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nRow > 0) {
    miisinfo.constrlb = bxCreateLogicalMatrix(nRow, 1);
    miisinfo.construb = bxCreateLogicalMatrix(nRow, 1);
    if (!miisinfo.constrlb || !miisinfo.construb) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nSos > 0) {
    miisinfo.sos = bxCreateLogicalMatrix(nSos, 1);
    if (!miisinfo.sos) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nIndicator > 0) {
    miisinfo.indicator = bxCreateLogicalMatrix(nIndicator, 1);
    if (!miisinfo.indicator) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nCol > 0) {
    ciisinfo.colLowerIIS = (int *) bxCalloc(nCol, sizeof(int));
    ciisinfo.colUpperIIS = (int *) bxCalloc(nCol, sizeof(int));
    if (!ciisinfo.colLowerIIS || !ciisinfo.colUpperIIS) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nRow > 0) {
    ciisinfo.rowLowerIIS = (int *) bxCalloc(nRow, sizeof(int));
    ciisinfo.rowUpperIIS = (int *) bxCalloc(nRow, sizeof(int));
    if (!ciisinfo.rowLowerIIS || !ciisinfo.rowUpperIIS) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nSos > 0) {
    ciisinfo.sosIIS = (int *) bxCalloc(nSos, sizeof(int));
    if (!ciisinfo.sosIIS) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nIndicator > 0) {
    ciisinfo.indicatorIIS = (int *) bxCalloc(nIndicator, sizeof(int));
    if (!ciisinfo.indicatorIIS) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ISMINIIS, &isMinIIS));

  if (nCol > 0) {
    COPTMEX_CALL(COPT_GetColLowerIIS(prob, nCol, NULL, ciisinfo.colLowerIIS));
    COPTMEX_CALL(COPT_GetColUpperIIS(prob, nCol, NULL, ciisinfo.colUpperIIS));
  }

  if (nRow > 0) {
    COPTMEX_CALL(COPT_GetRowLowerIIS(prob, nRow, NULL, ciisinfo.rowLowerIIS));
    COPTMEX_CALL(COPT_GetRowUpperIIS(prob, nRow, NULL, ciisinfo.rowUpperIIS));
  }

  if (nSos > 0) {
    COPTMEX_CALL(COPT_GetSOSIIS(prob, nSos, NULL, ciisinfo.sosIIS));
  }

  if (nIndicator > 0) {
    COPTMEX_CALL(COPT_GetIndicatorIIS(prob, nIndicator, NULL, ciisinfo.indicatorIIS));
  }

  *bxGetDoubles(miisinfo.isminiis) = isMinIIS;

  if (nCol > 0) {
    bxLogical *colLowerIIS_data = bxGetLogicals(miisinfo.varlb);
    bxLogical *colUpperIIS_data = bxGetLogicals(miisinfo.varub);
    for (int i = 0; i < nCol; ++i) {
      colLowerIIS_data[i] = ciisinfo.colLowerIIS[i];
      colUpperIIS_data[i] = ciisinfo.colUpperIIS[i];
    }
    bxFree(ciisinfo.colLowerIIS);
    bxFree(ciisinfo.colUpperIIS);
  }

  if (nRow > 0) {
    bxLogical *rowLowerIIS_data = bxGetLogicals(miisinfo.constrlb);
    bxLogical *rowUpperIIS_data = bxGetLogicals(miisinfo.construb);
    for (int i = 0; i < nRow; ++i) {
      rowLowerIIS_data[i] = ciisinfo.rowLowerIIS[i];
      rowUpperIIS_data[i] = ciisinfo.rowUpperIIS[i];
    }
    bxFree(ciisinfo.rowLowerIIS);
    bxFree(ciisinfo.rowUpperIIS);
  }

  if (nSos > 0) {
    bxLogical *sosIIS_data = bxGetLogicals(miisinfo.sos);
    for (int i = 0; i < nSos; ++i) {
      sosIIS_data[i] = ciisinfo.sosIIS[i];
    }
    bxFree(ciisinfo.sosIIS);
  }

  if (nIndicator > 0) {
    bxLogical *indicatorIIS_data = bxGetLogicals(miisinfo.indicator);
    for (int i = 0; i < nIndicator; ++i) {
      indicatorIIS_data[i] = ciisinfo.indicatorIIS[i];
    }
    bxFree(ciisinfo.indicatorIIS);
  }

  iisInfo = bxCreateStructMatrix(1, 1, 0, NULL);
  if (!iisInfo) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // 'isminiis'
  bxAddField(iisInfo, COPTMEX_IIS_ISMINIIS);
  bxSetField(iisInfo, 0, COPTMEX_IIS_ISMINIIS, miisinfo.isminiis);

  if (nCol > 0) {
    // 'varlb'
    bxAddField(iisInfo, COPTMEX_IIS_VARLB);
    bxSetField(iisInfo, 0, COPTMEX_IIS_VARLB, miisinfo.varlb);
    // 'varub'
    bxAddField(iisInfo, COPTMEX_IIS_VARUB);
    bxSetField(iisInfo, 0, COPTMEX_IIS_VARUB, miisinfo.varub);
  }

  if (nRow > 0) {
    // 'constrlb'
    bxAddField(iisInfo, COPTMEX_IIS_CONSTRLB);
    bxSetField(iisInfo, 0, COPTMEX_IIS_CONSTRLB, miisinfo.constrlb);
    // 'construb'
    bxAddField(iisInfo, COPTMEX_IIS_CONSTRUB);
    bxSetField(iisInfo, 0, COPTMEX_IIS_CONSTRUB, miisinfo.construb);
  }

  if (nSos > 0) {
    // 'sos'
    bxAddField(iisInfo, COPTMEX_IIS_SOS);
    bxSetField(iisInfo, 0, COPTMEX_IIS_SOS, miisinfo.sos);
  }

  if (nIndicator > 0) {
    // 'indicator'
    bxAddField(iisInfo, COPTMEX_IIS_INDICATOR);
    bxSetField(iisInfo, 0, COPTMEX_IIS_INDICATOR, miisinfo.indicator);
  }

  // Write out IIS problem
  COPTMEX_CALL(COPT_WriteIIS(prob, "result.iis"));

  *out_iis = iisInfo;

exit_cleanup:
  if (retcode != COPT_RETCODE_OK) {
    *out_iis = NULL;
  }

  return retcode;
}

/* Compute IIS for infeasible problem */
int COPTMEX_computeIIS(copt_prob *prob, bxArray **out_iis, int ifRetResult) {
  int retcode = COPT_RETCODE_OK;
  int modelStatus = 0;
  int isMIP = 0;

  // Try to find IIS for the given problem
  COPTMEX_CALL(COPT_ComputeIIS(prob));

  // Extract IIS information
  if (ifRetResult == 1) {
    COPTMEX_CALL(COPTMEX_getIIS(prob, out_iis));
  }

exit_cleanup:
  return retcode;
}

/* Extract feasibility relaxation information */
static int COPTMEX_getFeasRelax(copt_prob *prob, bxArray **out_relax) {
  int retcode = COPT_RETCODE_OK;
  int nRow = 0, nCol = 0;
  int hasFeasRelax = 0;
  bxArray *relaxInfo = NULL;
  coptmex_crelaxinfo crelaxinfo;
  coptmex_mrelaxinfo mrelaxinfo;

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_HASFEASRELAXSOL, &hasFeasRelax));
  if (hasFeasRelax == 0) {
    *out_relax = NULL;
    goto exit_cleanup;
  }

  COPTMEX_initCRelaxInfo(&crelaxinfo);
  COPTMEX_initMRelaxInfo(&mrelaxinfo);

  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &nRow));
  COPTMEX_CALL(COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &nCol));

  mrelaxinfo.relaxobj = bxCreateDoubleMatrix(1, 1, bxREAL);
  if (!mrelaxinfo.relaxobj) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  if (nCol > 0) {
    mrelaxinfo.relaxlb = bxCreateDoubleMatrix(nCol, 1, bxREAL);
    mrelaxinfo.relaxub = bxCreateDoubleMatrix(nCol, 1, bxREAL);
    if (!mrelaxinfo.relaxlb || !mrelaxinfo.relaxub) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nRow > 0) {
    mrelaxinfo.relaxlhs = bxCreateDoubleMatrix(nRow, 1, bxREAL);
    mrelaxinfo.relaxrhs = bxCreateDoubleMatrix(nRow, 1, bxREAL);
    if (!mrelaxinfo.relaxlhs || !mrelaxinfo.relaxrhs) {
      retcode = COPT_RETCODE_MEMORY;
      goto exit_cleanup;
    }
  }

  if (nCol > 0) {
    crelaxinfo.colLowRlx = bxGetDoubles(mrelaxinfo.relaxlb);
    crelaxinfo.colUppRlx = bxGetDoubles(mrelaxinfo.relaxub);
  }

  if (nRow > 0) {
    crelaxinfo.rowLowRlx = bxGetDoubles(mrelaxinfo.relaxlhs);
    crelaxinfo.rowUppRlx = bxGetDoubles(mrelaxinfo.relaxrhs);
  }

  COPTMEX_CALL(COPT_GetDblAttr(prob, COPT_DBLATTR_FEASRELAXOBJ, &crelaxinfo.dObjVal));
  *bxGetDoubles(mrelaxinfo.relaxobj) = crelaxinfo.dObjVal;

  if (nCol > 0) {
    COPTMEX_CALL(COPT_GetColInfo(prob, COPT_DBLINFO_RELAXLB, nCol, NULL, crelaxinfo.colLowRlx));
    COPTMEX_CALL(COPT_GetColInfo(prob, COPT_DBLINFO_RELAXUB, nCol, NULL, crelaxinfo.colUppRlx));
  }

  if (nRow > 0) {
    COPTMEX_CALL(COPT_GetRowInfo(prob, COPT_DBLINFO_RELAXLB, nRow, NULL, crelaxinfo.rowLowRlx));
    COPTMEX_CALL(COPT_GetRowInfo(prob, COPT_DBLINFO_RELAXUB, nRow, NULL, crelaxinfo.rowUppRlx));
  }

  relaxInfo = bxCreateStructMatrix(1, 1, 0, NULL);
  if (!relaxInfo) {
    retcode = COPT_RETCODE_MEMORY;
    goto exit_cleanup;
  }

  // 'relaxobj'
  bxAddField(relaxInfo, COPTMEX_FEASRELAX_OBJ);
  bxSetField(relaxInfo, 0, COPTMEX_FEASRELAX_OBJ, mrelaxinfo.relaxobj);

  if (nCol > 0) {
    // 'relaxlb'
    bxAddField(relaxInfo, COPTMEX_FEASRELAX_LB);
    bxSetField(relaxInfo, 0, COPTMEX_FEASRELAX_LB, mrelaxinfo.relaxlb);
    // 'relaxub'
    bxAddField(relaxInfo, COPTMEX_FEASRELAX_UB);
    bxSetField(relaxInfo, 0, COPTMEX_FEASRELAX_UB, mrelaxinfo.relaxub);
  }

  if (nRow > 0) {
    // 'relaxlhs'
    bxAddField(relaxInfo, COPTMEX_FEASRELAX_LHS);
    bxSetField(relaxInfo, 0, COPTMEX_FEASRELAX_LHS, mrelaxinfo.relaxlhs);
    // 'relaxrhs'
    bxAddField(relaxInfo, COPTMEX_FEASRELAX_RHS);
    bxSetField(relaxInfo, 0, COPTMEX_FEASRELAX_RHS, mrelaxinfo.relaxrhs);
  }

  // Write out feasibility relaxation problem
  COPTMEX_CALL(COPT_WriteRelax(prob, "result.relax"));

  *out_relax = relaxInfo;

exit_cleanup:
  if (retcode != COPT_RETCODE_OK) {
    *out_relax = NULL;
  }

  return retcode;
}

int COPTMEX_feasRelax(copt_prob *prob, const bxArray *penalty, bxArray **out_relax, int ifRetResult) {
  int retcode = COPT_RETCODE_OK;

  bxArray *lbpen = NULL;
  bxArray *ubpen = NULL;
  bxArray *rhspen = NULL;
  bxArray *upppen = NULL;

  double *colLowPen = NULL;
  double *colUppPen = NULL;
  double *rowBndPen = NULL;
  double *rowUppPen = NULL;

  if (COPTMEX_checkPenalty(prob, penalty) == 0) {
    goto exit_cleanup;
  }

  lbpen = bxGetField(penalty, 0, COPTMEX_PENALTY_LBPEN);
  ubpen = bxGetField(penalty, 0, COPTMEX_PENALTY_UBPEN);
  rhspen = bxGetField(penalty, 0, COPTMEX_PENALTY_RHSPEN);
  upppen = bxGetField(penalty, 0, COPTMEX_PENALTY_UPPPEN);

  if (lbpen != NULL) {
    colLowPen = bxGetDoubles(lbpen);
  }
  if (ubpen != NULL) {
    colUppPen = bxGetDoubles(ubpen);
  }
  if (rhspen != NULL) {
    rowBndPen = bxGetDoubles(rhspen);
  }
  if (upppen != NULL) {
    rowUppPen = bxGetDoubles(upppen);
  }

  // Compute the feasibility relaxation
  COPTMEX_CALL(COPT_FeasRelax(prob, colLowPen, colUppPen, rowBndPen, rowUppPen));

  // Extract feasibility relaxation information
  if (ifRetResult == 1) {
    COPTMEX_CALL(COPTMEX_getFeasRelax(prob, out_relax));
  }

exit_cleanup:
  return retcode;
}
