#include "coptmex.h"
#include <stdlib.h>

/* Initialize MEX-style version */
static void COPTMEX_initVersion(coptmex_mversion *version) {
  version->major     = NULL;
  version->minor     = NULL;
  version->technical = NULL;
  return;
}

/* Initialize C-style problem */
static void COPTMEX_initCProb(coptmex_cprob *cprob) {
  // The main part of problem
  cprob->nRow       = 0;
  cprob->nCol       = 0;
  cprob->nElem      = 0;
  cprob->nObjSen    = COPT_MINIMIZE;
  cprob->dObjConst  = 0.0;

  cprob->colMatBeg  = NULL;
  cprob->colMatIdx  = NULL;
  cprob->colMatElem = NULL;

  cprob->colCost    = NULL;
  cprob->colLower   = NULL;
  cprob->colUpper   = NULL;
  cprob->rowLower   = NULL;
  cprob->rowUpper   = NULL;

  cprob->colType    = NULL;
  cprob->rowSense   = NULL;
  cprob->colNames   = NULL;
  cprob->rowNames   = NULL;

  // The optional SOS part
  cprob->nSos       = 0;
  cprob->nSosSize   = 0;
  cprob->sosType    = NULL;
  cprob->sosMatBeg  = NULL;
  cprob->sosMatCnt  = NULL;
  cprob->sosMatIdx  = NULL;
  cprob->sosMatWt   = NULL;

  // The optional indicator part
  cprob->nIndicator = 0;

  // The optional cone part
  cprob->nCone      = 0;
  cprob->nConeSize  = 0;
  cprob->coneType   = NULL;
  cprob->coneBeg    = NULL;
  cprob->coneCnt    = NULL;
  cprob->coneIdx    = NULL;

  // The optional Q objective part
  cprob->nQElem     = 0;

  // The optional quadratic constraint part
  cprob->nQConstr   = 0;

  // The optional advanced information
  cprob->hasBasis   = 0;
  cprob->colBasis   = NULL;
  cprob->rowBasis   = NULL;
  return;
}

/* Initialize MEX-style problem */
static void COPTMEX_initMProb(coptmex_mprob *mprob) {
  // The main part of problem
  mprob->objsen      = NULL;
  mprob->objcon      = NULL;
  mprob->A           = NULL;
  mprob->obj         = NULL;
  mprob->lb          = NULL;
  mprob->ub          = NULL;
  mprob->vtype       = NULL;
  mprob->varnames    = NULL;
  mprob->sense       = NULL;
  mprob->lhs         = NULL;
  mprob->rhs         = NULL;
  mprob->constrnames = NULL;

  // The optional SOS part
  mprob->sos         = NULL;

  // The optional indicator part
  mprob->indicator   = NULL;

  // The optional cone part
  mprob->cone        = NULL;

  // The optional Q objective part
  mprob->qobj        = NULL;

  // The optional quadratic constraint part
  mprob->quadcon     = NULL;

  // The optional advanced information
  mprob->varbasis    = NULL;
  mprob->constrbasis = NULL;

  mprob->value       = NULL;
  mprob->slack       = NULL;
  mprob->dual        = NULL;
  mprob->redcost     = NULL;

  mprob->mipstart    = NULL;
  return;
}

/* Initialize C-style cone problem */
static void COPTMEX_initCConeProb(coptmex_cconeprob *cconeprob) {
  cconeprob->nCol          = 0;
  cconeprob->nRow          = 0;

  cconeprob->nObjSense     = COPT_MINIMIZE;
  cconeprob->dObjConst     = 0.0;

  cconeprob->nFree         = 0;
  cconeprob->nPositive     = 0;
  cconeprob->nCone         = 0;
  cconeprob->nRotateCone   = 0;
  cconeprob->nPSD          = 0;

  cconeprob->coneDim       = NULL;
  cconeprob->rotateConeDim = NULL;
  cconeprob->psdDim        = NULL;

  cconeprob->colObj        = NULL;

  cconeprob->nQObjElem     = 0;
  cconeprob->qObjRow       = NULL;
  cconeprob->qObjCol       = NULL;
  cconeprob->qObjElem      = NULL;

  cconeprob->colMatBeg     = NULL;
  cconeprob->colMatIdx     = NULL;
  cconeprob->colMatElem    = NULL;

  cconeprob->rowRhs        = NULL;
}

/* Initialize MEX-style cone problem */
static void COPTMEX_initMConeProb(coptmex_mconeprob *mconeprob) {
  mconeprob->c      = NULL;
  mconeprob->A      = NULL;
  mconeprob->b      = NULL;

  mconeprob->K      = NULL;
  mconeprob->f      = NULL;
  mconeprob->l      = NULL;
  mconeprob->q      = NULL;
  mconeprob->r      = NULL;
  mconeprob->s      = NULL;

  mconeprob->objsen = NULL;
  mconeprob->objcon = NULL;
  mconeprob->Q      = NULL;
}

/* Initialize C-style LP solution */
static void COPTMEX_initCLpSol(coptmex_clpsol *clpsol) {
  clpsol->nRow         = 0;
  clpsol->nCol         = 0;
  clpsol->nPSD         = 0;
  clpsol->nPSDLen      = 0;
  clpsol->nPSDConstr   = 0;
  clpsol->nQConstr     = 0;
  clpsol->hasBasis     = 0;
  clpsol->hasLpSol     = 0;

  clpsol->nStatus      = COPT_LPSTATUS_UNSTARTED;
  clpsol->nSimplexIter = 0;
  clpsol->nBarrierIter = 0;
  clpsol->dSolvingTime = 0.0;
  clpsol->dObjVal      = COPT_INFINITY;

  clpsol->colBasis     = NULL;
  clpsol->rowBasis     = NULL;
  clpsol->colValue     = NULL;
  clpsol->colDual      = NULL;
  clpsol->rowSlack     = NULL;
  clpsol->rowDual      = NULL;
  clpsol->primalRay    = NULL;
  clpsol->dualFarkas   = NULL;

  clpsol->qRowSlack    = NULL;

  clpsol->psdColValue  = NULL;
  clpsol->psdColDual   = NULL;
  clpsol->psdRowSlack  = NULL;
  clpsol->psdRowDual   = NULL;
  return;
}

/* Initialize C-style MIP solution */
static void COPTMEX_initCMipSol(coptmex_cmipsol *cmipsol) {
  cmipsol->nRow         = 0;
  cmipsol->nCol         = 0;
  cmipsol->hasMipSol    = 0;

  cmipsol->nStatus      = COPT_MIPSTATUS_UNSTARTED;
  cmipsol->nSimplexIter = 0;
  cmipsol->nNodeCnt     = 0;
  cmipsol->dBestGap     = COPT_INFINITY;
  cmipsol->dSolvingTime = 0.0;
  cmipsol->dObjVal      = COPT_INFINITY;
  cmipsol->dBestBnd     = -COPT_INFINITY;

  cmipsol->colValue     = NULL;

  cmipsol->nSolPool     = 0;
  return;
}

/* Initialize C-style IIS information */
static void COPTMEX_initCIISInfo(coptmex_ciisinfo *ciisinfo) {
  ciisinfo->isMinIIS     = 0;
  ciisinfo->colLowerIIS  = NULL;
  ciisinfo->colUpperIIS  = NULL;
  ciisinfo->rowLowerIIS  = NULL;
  ciisinfo->rowUpperIIS  = NULL;
  ciisinfo->sosIIS       = NULL;
  ciisinfo->indicatorIIS = NULL;
  return;
}

/* Initialize C-style feasibility relaxation information */
static void COPTMEX_initCRelaxInfo(coptmex_crelaxinfo *crelaxinfo) {
  crelaxinfo->dObjVal       = 0.0;
  crelaxinfo->colLowRlx = NULL;
  crelaxinfo->colUppRlx = NULL;
  crelaxinfo->rowLowRlx = NULL;
  crelaxinfo->rowUppRlx = NULL;
  return;
}

/* Initialize MEX-style LP solution */
static void COPTMEX_initMLpSol(coptmex_mlpsol *mlpsol) {
  mlpsol->status      = NULL;
  mlpsol->simplexiter = NULL;
  mlpsol->barrieriter = NULL;
  mlpsol->solvingtime = NULL;
  mlpsol->objval      = NULL;
  mlpsol->varbasis    = NULL;
  mlpsol->constrbasis = NULL;
  mlpsol->value       = NULL;
  mlpsol->redcost     = NULL;
  mlpsol->slack       = NULL;
  mlpsol->dual        = NULL;
  mlpsol->ray         = NULL;
  mlpsol->farkas      = NULL;
  mlpsol->qcslack     = NULL;
  mlpsol->psdcolvalue = NULL;
  mlpsol->psdcoldual  = NULL;
  mlpsol->psdrowslack = NULL;
  mlpsol->psdrowdual  = NULL;
  return;
}

/* Initialize MEX-style MIP solution */
static void COPTMEX_initMMipSol(coptmex_mmipsol *mmipsol) {
  mmipsol->status      = NULL;
  mmipsol->simplexiter = NULL;
  mmipsol->nodecnt     = NULL;
  mmipsol->bestgap     = NULL;
  mmipsol->solvingtime = NULL;
  mmipsol->objval      = NULL;
  mmipsol->bestbnd     = NULL;
  mmipsol->value       = NULL;

  mmipsol->solpool     = NULL;
  return;
}

/* Initialize MEX-style IIS information */
static void COPTMEX_initMIISInfo(coptmex_miisinfo *miisinfo) {
  miisinfo->isminiis  = NULL;
  miisinfo->varlb     = NULL;
  miisinfo->varub     = NULL;
  miisinfo->constrlb  = NULL;
  miisinfo->construb  = NULL;
  miisinfo->sos       = NULL;
  miisinfo->indicator = NULL;
  return;
}

/* Initialize MEX-style feasibility relaxation information */
static void COPTMEX_initMRelaxInfo(coptmex_mrelaxinfo *mrelaxinfo) {
  mrelaxinfo->relaxobj = NULL;
  mrelaxinfo->relaxlb  = NULL;
  mrelaxinfo->relaxub  = NULL;
  mrelaxinfo->relaxlhs = NULL;
  mrelaxinfo->relaxrhs = NULL;
  return;
}

/* Check parts of penalty */
static int COPTMEX_checkPenalty(copt_prob *prob, const bxArray *penalty) {
  int isvalid = 1;
  char msgbuf[COPT_BUFFSIZE];

  bxArray *lbpen = NULL;
  bxArray *ubpen = NULL;
  bxArray *rhspen = NULL;
  bxArray *upppen = NULL;

  int nCol = 0, nRow = 0;

  COPT_GetIntAttr(prob, COPT_INTATTR_COLS, &nCol);
  COPT_GetIntAttr(prob, COPT_INTATTR_ROWS, &nRow);

  lbpen = bxGetField(penalty, 0, COPTMEX_PENALTY_LBPEN);
  ubpen = bxGetField(penalty, 0, COPTMEX_PENALTY_UBPEN);
  rhspen = bxGetField(penalty, 0, COPTMEX_PENALTY_RHSPEN);
  upppen = bxGetField(penalty, 0, COPTMEX_PENALTY_UPPPEN);

  if (lbpen != NULL) {
    if (!bxIsDouble(lbpen)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_LBPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(lbpen) != nCol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_LBPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }

  if (ubpen != NULL) {
    if (!bxIsDouble(ubpen)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_UBPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(ubpen) != nCol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_UBPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }

  if (rhspen != NULL) {
    if (!bxIsDouble(rhspen)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_RHSPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(rhspen) != nRow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_RHSPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }

  if (upppen != NULL) {
    if (!bxIsDouble(upppen)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_UPPPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(upppen) != nRow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "penalty.%s", COPTMEX_PENALTY_UPPPEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }

exit_cleanup:
  return isvalid;
}

/* Check all parts of a cone problem */
static int COPTMEX_checkConeModel(bxArray *conedata) {
  int nrow = 0, ncol = 0;
  int isvalid = 1;
  char msgbuf[COPT_BUFFSIZE];

  // 'conedata'
  if (conedata != NULL) {
    if (!bxIsStruct(conedata)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_CONEDATA);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }

    // 'A'
    bxArray *A = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_A);
    if (A == NULL) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_A);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
      goto exit_cleanup;
    } else {
      if (!bxIsSparse(A)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_A);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      } else {
        nrow = (int) bxGetM(A);
        ncol = (int) bxGetN(A);
      }
    }

    // 'K'
    bxArray *K = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_K);
    if (K == NULL) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_K);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
      goto exit_cleanup;
    } else {
      if (!bxIsStruct(K)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_K);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 'f'
    bxArray *f = bxGetField(conedata, 0, COPTMEX_MODEL_CONEK_F);
    if (f != NULL) {
      if (!bxIsScalar(f) || bxIsChar(f)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONEK_F);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 'l'
    bxArray *l = bxGetField(conedata, 0, COPTMEX_MODEL_CONEK_L);
    if (l != NULL) {
      if (!bxIsScalar(l) || bxIsChar(l)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONEK_L);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 'q'
    bxArray *q = bxGetField(conedata, 0, COPTMEX_MODEL_CONEK_Q);
    if (q != NULL) {
      if (!bxIsDouble(q)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONEK_Q);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 'r'
    bxArray *r = bxGetField(conedata, 0, COPTMEX_MODEL_CONEK_R);
    if (r != NULL) {
      if (!bxIsDouble(r)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONEK_R);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 's'
    bxArray *s = bxGetField(conedata, 0, COPTMEX_MODEL_CONEK_S);
    if (s != NULL) {
      if (!bxIsDouble(s)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONEK_S);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 'c'
    bxArray *c = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_C);
    if (c != NULL) {
      if (!bxIsDouble(c)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_C);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      } else {
        if (bxGetNumberOfElements(c) != ncol) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_C);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
          goto exit_cleanup;
        }
      }
    }

    // 'b'
    bxArray *b = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_B);
    if (b != NULL) {
      if (!bxIsDouble(b)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_B);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      } else {
        if (bxGetNumberOfElements(b) != nrow) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.conedata.%s", COPTMEX_MODEL_CONE_B);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
          goto exit_cleanup;
        }
      }
    }

    // 'objsen'
    bxArray *objsen = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_OBJSEN);
    if (objsen != NULL) {
      if (!bxIsChar(objsen)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_CONE_OBJSEN);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }

    // 'objcon'
    bxArray *objcon = bxGetField(conedata, 0, COPTMEX_MODEL_CONE_OBJCON);
    if (objcon != NULL) {
      if (!bxIsScalar(objcon) || bxIsChar(objcon)) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_CONE_OBJCON);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
        goto exit_cleanup;
      }
    }
  }

exit_cleanup:
  return isvalid;
}

/* Check all parts of a problem */
static int COPTMEX_checkModel(coptmex_mprob *mprob) {
  int nrow = 0, ncol = 0;
  int isvalid = 1;
  char msgbuf[COPT_BUFFSIZE];

  // 'objsen'
  if (mprob->objsen != NULL) {
    if (!bxIsChar(mprob->objsen)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_OBJSEN);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'objcon'
  if (mprob->objcon != NULL) {
    if (!bxIsScalar(mprob->objcon) || bxIsChar(mprob->objcon)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_OBJCON);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'A'
  if (mprob->A == NULL) {
    isvalid = 0;
    snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_A);
    COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
    goto exit_cleanup;
  } else {
    if (!bxIsSparse(mprob->A)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_A);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    } else {
      nrow = (int) bxGetM(mprob->A);
      ncol = (int) bxGetN(mprob->A);
    }
  }
  // 'obj'
  if (mprob->obj != NULL) {
    if (!bxIsDouble(mprob->obj)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_OBJ);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->obj) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_OBJ);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'lb'
  if (mprob->lb != NULL) {
    if (!bxIsDouble(mprob->lb)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_LB);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->lb) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_LB);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'ub'
  if (mprob->ub != NULL) {
    if (!bxIsDouble(mprob->ub)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_UB);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->ub) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_UB);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'vtype'
  if (mprob->vtype != NULL) {
    if (!bxIsChar(mprob->vtype)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_VTYPE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->vtype) != ncol &&
        bxGetNumberOfElements(mprob->vtype) != 1) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_VTYPE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'sense'
  if (mprob->sense == NULL) {
    // 'lhs'
    if (!bxIsDouble(mprob->lhs)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_LHS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->lhs) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_LHS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }

    // 'rhs'
    if (!bxIsDouble(mprob->rhs)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_RHS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->rhs) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_RHS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  } else {
    // 'sense'
    if (!bxIsChar(mprob->sense)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_SENSE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->sense) != nrow &&
        bxGetNumberOfElements(mprob->sense) != 1) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_SENSE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }

    // 'rhs'
    if (!bxIsDouble(mprob->rhs)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_RHS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->rhs) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_RHS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'varnames'
  if (mprob->varnames != NULL) {
    if (!bxIsCell(mprob->varnames)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_VARNAME);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->varnames) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_VARNAME);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }

    for (int i = 0; i < bxGetNumberOfElements(mprob->varnames); ++i) {
      if (!bxIsChar(bxGetCell(mprob->varnames, i))) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s{%d}", COPTMEX_MODEL_VARNAME, i);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      }
    }
  }
  // 'constrnames'
  if (mprob->constrnames != NULL) {
    if (!bxIsCell(mprob->constrnames)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_CONNAME);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->constrnames) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_CONNAME);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }

    for (int i = 0; i < bxGetNumberOfElements(mprob->constrnames); ++i) {
      if (!bxIsChar(bxGetCell(mprob->constrnames, i))) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s{%d}", COPTMEX_MODEL_CONNAME, i);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      }
    }
  }

  // 'sos'
  if (mprob->sos != NULL) {
    if (!bxIsStruct(mprob->sos)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_SOS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }

    for (int i = 0; i < bxGetNumberOfElements(mprob->sos); ++i) {
      bxArray *sostype = bxGetField(mprob->sos, i, COPTMEX_MODEL_SOSTYPE);
      bxArray *sosvars = bxGetField(mprob->sos, i, COPTMEX_MODEL_SOSVARS);
      bxArray *soswght = bxGetField(mprob->sos, i, COPTMEX_MODEL_SOSWEIGHT);

      if (sostype == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_SOS,
                 i, COPTMEX_MODEL_SOSTYPE);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsScalar(sostype) || bxIsChar(sostype)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_SOS,
                   i, COPTMEX_MODEL_SOSTYPE);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (sosvars == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_SOS,
                 i, COPTMEX_MODEL_SOSVARS);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsDouble(sosvars) || bxIsScalar(sosvars) || bxIsSparse(sosvars)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_SOS,
                   i, COPTMEX_MODEL_SOSVARS);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (soswght != NULL) {
        if (!bxIsDouble(soswght) || bxIsScalar(soswght) || bxIsSparse(soswght)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_SOS,
                   i, COPTMEX_MODEL_SOSWEIGHT);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }

        if (bxGetNumberOfElements(sosvars) != bxGetNumberOfElements(soswght)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s and problem.%s(%d).%s",
                   COPTMEX_MODEL_SOS, i, COPTMEX_MODEL_SOSVARS,
                   COPTMEX_MODEL_SOS, i, COPTMEX_MODEL_SOSWEIGHT);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
          goto exit_cleanup;
        }
      }
    }
  }

  // 'indicator'
  if (mprob->indicator != NULL) {
    if (!bxIsStruct(mprob->indicator)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_INDICATOR);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }

    for (int i = 0; i < bxGetNumberOfElements(mprob->indicator); ++i) {
      bxArray *binvar = bxGetField(mprob->indicator, i, COPTMEX_MODEL_INDICBINVAR);
      bxArray *binval = bxGetField(mprob->indicator, i, COPTMEX_MODEL_INDICBINVAL);
      bxArray *indicA = bxGetField(mprob->indicator, i, COPTMEX_MODEL_INDICROW);
      bxArray *rSense = bxGetField(mprob->indicator, i, COPTMEX_MODEL_INDICSENSE);
      bxArray *rowBnd = bxGetField(mprob->indicator, i, COPTMEX_MODEL_INDICRHS);

      if (binvar == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                 i, COPTMEX_MODEL_INDICBINVAR);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsScalar(binvar) || bxIsChar(binvar)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                   i, COPTMEX_MODEL_INDICBINVAR);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (binval == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                 i, COPTMEX_MODEL_INDICBINVAL);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsScalar(binval) || bxIsChar(binval)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                   i, COPTMEX_MODEL_INDICBINVAL);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (indicA == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                 i, COPTMEX_MODEL_INDICROW);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsDouble(indicA)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                   i, COPTMEX_MODEL_INDICROW);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
        if (!bxIsSparse(indicA)) {
          if (bxGetNumberOfElements(indicA) != ncol) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                     i, COPTMEX_MODEL_INDICROW);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
            goto exit_cleanup;
          }
        } else {
          if (bxGetN(indicA) != 1) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                     i, COPTMEX_MODEL_INDICROW);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
            goto exit_cleanup;
          }
        }
      }

      if (rSense == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                 i, COPTMEX_MODEL_INDICSENSE);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsChar(rSense) || !bxIsScalar(rSense)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                   i, COPTMEX_MODEL_INDICSENSE);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (rowBnd == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                 i, COPTMEX_MODEL_INDICRHS);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsScalar(rowBnd) || bxIsChar(rowBnd)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_INDICATOR,
                   i, COPTMEX_MODEL_INDICRHS);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }
    }
  }

  // 'cone'
  if (mprob->cone != NULL) {
    if (!bxIsStruct(mprob->cone)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_CONE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }

    for (int i = 0; i < bxGetNumberOfElements(mprob->cone); ++i) {
      bxArray *conetype = bxGetField(mprob->cone, i, COPTMEX_MODEL_CONETYPE);
      bxArray *conevars = bxGetField(mprob->cone, i, COPTMEX_MODEL_CONEVARS);

      if (conetype == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_CONE,
                 i, COPTMEX_MODEL_CONETYPE);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsScalar(conetype) || bxIsChar(conetype)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_CONE,
                   i, COPTMEX_MODEL_CONETYPE);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (conevars == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_CONE,
                 i, COPTMEX_MODEL_CONEVARS);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (!bxIsDouble(conevars) || bxIsScalar(conevars) || bxIsSparse(conevars)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_CONE,
                   i, COPTMEX_MODEL_CONEVARS);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }
    }
  }

  // 'Q'
  if (mprob->qobj != NULL) {
    if (!bxIsSparse(mprob->qobj)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_QUADOBJ);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetM(mprob->qobj) != ncol || bxGetN(mprob->qobj) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_QUADOBJ);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }

  // 'quadcon'
  if (mprob->quadcon != NULL) {
    if (!bxIsStruct(mprob->quadcon)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_MODEL_QUADCON);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }

    for (int i = 0; i < bxGetNumberOfElements(mprob->quadcon); ++i) {
      bxArray *QcMat = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCSPMAT);
      bxArray *QcRow = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCROW);
      bxArray *QcCol = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCCOL);
      bxArray *QcVal = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCVAL);
      bxArray *QcLinear = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCLINEAR);
      bxArray *QcSense = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCSENSE);
      bxArray *QcRhs = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCRHS);
      bxArray *QcName = bxGetField(mprob->quadcon, i, COPTMEX_MODEL_QCNAME);

      if (QcMat != NULL) {
        if (!bxIsSparse(QcMat)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                   i, COPTMEX_MODEL_QCSPMAT);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
          goto exit_cleanup;
        }
        if (bxGetM(QcMat) != bxGetN(QcMat)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                   i, COPTMEX_MODEL_QCSPMAT);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
          goto exit_cleanup;
        }
      } else {
        if (QcRow != NULL && QcCol != NULL && QcVal != NULL) {
          if (!bxIsDouble(QcRow)) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                     i, COPTMEX_MODEL_QCROW);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
            goto exit_cleanup;
          }
          if (!bxIsDouble(QcCol)) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                     i, COPTMEX_MODEL_QCCOL);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
            goto exit_cleanup;
          }
          if (!bxIsDouble(QcVal)) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                     i, COPTMEX_MODEL_QCVAL);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
            goto exit_cleanup;
          }
          if (bxGetNumberOfElements(QcRow) != bxGetNumberOfElements(QcCol)) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                     i, COPTMEX_MODEL_QCROW);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                     i, COPTMEX_MODEL_QCCOL);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
            goto exit_cleanup;
          }
          if (bxGetNumberOfElements(QcRow) != bxGetNumberOfElements(QcVal)) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                     i, COPTMEX_MODEL_QCVAL);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
            goto exit_cleanup;
          }
        }
      }

      if (QcLinear != NULL) {
        if (!bxIsDouble(QcLinear)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                   i, COPTMEX_MODEL_QCLINEAR);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
        if (!bxIsSparse(QcLinear)) {
          if (bxGetNumberOfElements(QcLinear) != ncol) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON,
                     i, COPTMEX_MODEL_QCLINEAR);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
            goto exit_cleanup;
          }
        } else {
          if (bxGetN(QcLinear) != 1) {
            isvalid = 0;
            snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON,
                     i, COPTMEX_MODEL_QCLINEAR);
            COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
            goto exit_cleanup;
          }
        }
      }

      if (QcMat == NULL && QcRow == NULL && QcCol == NULL && QcVal == NULL && QcLinear == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d)", COPTMEX_MODEL_QUADCON, i);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      }

      if (QcSense != NULL) {
        if (!bxIsChar(QcSense) || !bxIsScalar(QcSense)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                   i, COPTMEX_MODEL_QCSENSE);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (QcRhs == NULL) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                 i, COPTMEX_MODEL_QCRHS);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      } else {
        if (bxIsChar(QcRhs) || !bxIsScalar(QcRhs)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                   i, COPTMEX_MODEL_QCRHS);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }

      if (QcName != NULL) {
        if (!bxIsChar(QcName)) {
          isvalid = 0;
          snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s(%d).%s", COPTMEX_MODEL_QUADCON, 
                   i, COPTMEX_MODEL_QCNAME);
          COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
          goto exit_cleanup;
        }
      }
    }
  }

  // 'varbasis'
  if (mprob->varbasis != NULL) {
    if (!bxIsDouble(mprob->varbasis) || bxIsScalar(mprob->varbasis) ||
        bxIsSparse(mprob->varbasis)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_VARBASIS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->varbasis) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_VARBASIS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'constrbasis'
  if (mprob->constrbasis != NULL) {
    if (!bxIsDouble(mprob->constrbasis) || bxIsScalar(mprob->constrbasis) ||
        bxIsSparse(mprob->constrbasis)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_CONBASIS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->constrbasis) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_CONBASIS);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'x'
  if (mprob->value != NULL) {
    if (!bxIsDouble(mprob->value) || bxIsScalar(mprob->value) || bxIsSparse(mprob->value)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_VALUE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->value) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_VALUE);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'rc'
  if (mprob->redcost != NULL) {
    if (!bxIsDouble(mprob->redcost) || bxIsScalar(mprob->redcost) || bxIsSparse(mprob->redcost)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_REDCOST);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->redcost) != ncol) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_REDCOST);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'slack'
  if (mprob->slack != NULL) {
    if (!bxIsDouble(mprob->slack) || bxIsScalar(mprob->slack) || bxIsSparse(mprob->slack)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_SLACK);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->slack) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_SLACK);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'pi'
  if (mprob->dual != NULL) {
    if (!bxIsDouble(mprob->dual) || bxIsScalar(mprob->dual) || bxIsSparse(mprob->dual)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_DUAL);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (bxGetNumberOfElements(mprob->dual) != nrow) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_RESULT_DUAL);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
      goto exit_cleanup;
    }
  }
  // 'start'
  if (mprob->mipstart != NULL) {
    if (!bxIsDouble(mprob->mipstart)) {
      isvalid = 0;
      snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_ADVINFO_MIPSTART);
      COPTMEX_errorMsg(COPTMEX_ERROR_BAD_TYPE, msgbuf);
      goto exit_cleanup;
    }
    if (!bxIsSparse(mprob->mipstart)) {
      if (bxGetNumberOfElements(mprob->mipstart) != ncol) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_ADVINFO_MIPSTART);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_NUM, msgbuf);
        goto exit_cleanup;
      }
    } else {
      if (bxGetN(mprob->mipstart) != 1) {
        isvalid = 0;
        snprintf(msgbuf, COPT_BUFFSIZE, "problem.%s", COPTMEX_ADVINFO_MIPSTART);
        COPTMEX_errorMsg(COPTMEX_ERROR_BAD_DATA, msgbuf);
        goto exit_cleanup;
      }
    }
  }

exit_cleanup:
  return isvalid;
}
