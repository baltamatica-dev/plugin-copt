/* SPDX-License-Identifier:  */
#include "coptmex.h"
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>


/* workaround */

void mexErrMsgIdAndTxt(const char *errorid, const char *errormsg, ...) {
    va_list ap;
    
    va_start(ap, errormsg);
    bxPrintf(errormsg, ap);
    va_end(ap);
    
    bxErrMsgTxt(errorid);
}

extern const char* INTPARAM_field_names[];
const char* bxGetFieldNameByNumber(const bxArray *pm, int fieldnumber) {
    return INTPARAM_field_names[fieldnumber];
}

/**
 * @brief 返回数组中的元素个数.
 * 
 * 如数组维度为 2x5x5，则返回 30
 */
int bxGetNumberOfElements(const bxArray *ba) {
    // need mxGetDimensions, mxGetNumberOfDimensions
    return bxGetM(ba) * bxGetN(ba);
}

/**
 * @brief 为结构数组添加一个新成员名.
 * 
 */
int bxAddField(bxArray *pm, const char *fieldname) {
    return -1; // error
};


bool bxIsScalar(const bxArray* ba) {
    return 1 == (bxGetM(ba)*bxGetN(ba));
}

double bxGetScalar(const bxArray* ba) {
    return *bxGetDoubles(ba);
}

int mexEvalString(const char *command) {
    return 1; // error
};


const char* copt_version_help = R"(
copt 绑定

    copt 插件
    Github:
    LICENSE:    
    Copyright (c) 2022 Chengyu HAN

版本信息
    bex SDK: 2.2.1

)"; /* copt_version_help */
BALTAM_PLUGIN_FCN(copt_version) {
    bxPrintf(copt_version_help);
} /* copt_version */


/**
 * @brief [可选] 插件初始化函数.
 *
 * bxPluginInit 由 load_plugin(name, args...) 调用
 * 用于进行一些初始化工作
 *
 * @param nInit
 * @param pInit[]
 */
int bxPluginInit(int nInit, const bxArray* pInit[]) {
    return 0;
} /* bxPluginInit */

/**
 * @brief [可选] 插件终止时清理函数.
 *
 * bxPluginFini 由 unload_plugin() 调用
 * 用于进行一些清理工作
 */
int bxPluginFini() {
    return 0;
} /* bxPluginFini */

/**
 * @brief 【必选】 列出插件提供的函数.
 *
 * bxPluginFunctions 返回 指向函数列表的指针.
 */
bexfun_info_t * bxPluginFunctions() {
    // 已定义的插件函数个数
    constexpr size_t TOTAL_PLUGIN_FUNCTIONS = 6;
    bexfun_info_t* func_list_dyn = new bexfun_info_t[TOTAL_PLUGIN_FUNCTIONS + 1];

    size_t i = 0;
    func_list_dyn[i].name = "copt_version";
    func_list_dyn[i].ptr  = copt_version;
    func_list_dyn[i].help = copt_version_help;

    i++;
    func_list_dyn[i].name = "__copt_computeiis_impl";
    func_list_dyn[i].ptr  = copt_computeiis;
    func_list_dyn[i].help = nullptr;
    
    i++;
    func_list_dyn[i].name = "__copt_feasrelax_impl";
    func_list_dyn[i].ptr  = copt_feasrelax;
    func_list_dyn[i].help = nullptr;
    
    i++;
    func_list_dyn[i].name = "__copt_read_impl";
    func_list_dyn[i].ptr  = copt_read;
    func_list_dyn[i].help = nullptr;
    
    i++;
    func_list_dyn[i].name = "__copt_solve_impl";
    func_list_dyn[i].ptr  = copt_solve;
    func_list_dyn[i].help = nullptr;
    
    i++;
    func_list_dyn[i].name = "__copt_write_impl";
    func_list_dyn[i].ptr  = copt_write;
    func_list_dyn[i].help = nullptr;

    // 最后一个元素, `name` 字段必须为空字符串 `""`
    i++;
    func_list_dyn[i].name = "";
    func_list_dyn[i].ptr  = nullptr;
    func_list_dyn[i].help = nullptr;

    assert((TOTAL_PLUGIN_FUNCTIONS == i));
    return func_list_dyn;
} /* bxPluginFunctions */



const char* INTPARAM_field_names[] = {
    COPT_INTPARAM_LOGGING,
    COPT_INTPARAM_LOGTOCONSOLE,
    COPT_INTPARAM_PRESOLVE,
    COPT_INTPARAM_SCALING,
    COPT_INTPARAM_DUALIZE,
    COPT_INTPARAM_LPMETHOD,
    COPT_INTPARAM_REQFARKASRAY,
    COPT_INTPARAM_DUALPRICE,
    COPT_INTPARAM_DUALPERTURB,
    COPT_INTPARAM_CUTLEVEL,
    COPT_INTPARAM_ROOTCUTLEVEL,
    COPT_INTPARAM_TREECUTLEVEL,
    COPT_INTPARAM_ROOTCUTROUNDS,
    COPT_INTPARAM_NODECUTROUNDS,
    COPT_INTPARAM_HEURLEVEL,
    COPT_INTPARAM_ROUNDINGHEURLEVEL,
    COPT_INTPARAM_DIVINGHEURLEVEL,
    COPT_INTPARAM_SUBMIPHEURLEVEL,
    COPT_INTPARAM_STRONGBRANCHING,
    COPT_INTPARAM_CONFLICTANALYSIS,
    COPT_INTPARAM_NODELIMIT,
    COPT_INTPARAM_MIPTASKS,
    COPT_INTPARAM_BARHOMOGENEOUS,
    COPT_INTPARAM_BARORDER,
    COPT_INTPARAM_BARITERLIMIT,
    COPT_INTPARAM_THREADS,
    COPT_INTPARAM_BARTHREADS,
    COPT_INTPARAM_SIMPLEXTHREADS,
    COPT_INTPARAM_CROSSOVERTHREADS,
    COPT_INTPARAM_CROSSOVER,
    COPT_INTPARAM_SDPMETHOD,
    COPT_INTPARAM_IISMETHOD,
    COPT_INTPARAM_FEASRELAXMODE,
    COPT_INTPARAM_MIPSTARTMODE,
    COPT_INTPARAM_MIPSTARTNODELIMIT
};


const char* model_fields[] = {
  COPTMEX_MODEL_OBJSEN,
  COPTMEX_MODEL_OBJCON,
  COPTMEX_MODEL_A,
  COPTMEX_MODEL_OBJ,
  COPTMEX_MODEL_LB,
  COPTMEX_MODEL_UB,
  COPTMEX_MODEL_VTYPE,
  COPTMEX_MODEL_VARNAME,
  COPTMEX_MODEL_SENSE,
  COPTMEX_MODEL_LHS,
  COPTMEX_MODEL_RHS,
  COPTMEX_MODEL_CONNAME,

/* optional part */
  COPTMEX_MODEL_SOS,
  COPTMEX_MODEL_SOSTYPE,
  COPTMEX_MODEL_SOSVARS,
  COPTMEX_MODEL_SOSWEIGHT,

  COPTMEX_MODEL_INDICATOR,
  COPTMEX_MODEL_INDICBINVAR,
  COPTMEX_MODEL_INDICBINVAL,
  COPTMEX_MODEL_INDICROW,
  COPTMEX_MODEL_INDICSENSE,
  COPTMEX_MODEL_INDICRHS,

  COPTMEX_MODEL_QUADOBJ,

  COPTMEX_MODEL_QUADCON,
  COPTMEX_MODEL_QCSPMAT,
  COPTMEX_MODEL_QCROW,
  COPTMEX_MODEL_QCCOL,
  COPTMEX_MODEL_QCVAL,
  COPTMEX_MODEL_QCLINEAR,
  COPTMEX_MODEL_QCSENSE,
  COPTMEX_MODEL_QCRHS,
  COPTMEX_MODEL_QCNAME,

  COPTMEX_MODEL_CONE,
  COPTMEX_MODEL_CONETYPE,
  COPTMEX_MODEL_CONEVARS,

  COPTMEX_MODEL_CONEDATA,
  COPTMEX_MODEL_CONE_OBJSEN,
  COPTMEX_MODEL_CONE_OBJCON,
  COPTMEX_MODEL_CONE_C,
  COPTMEX_MODEL_CONE_A,
  COPTMEX_MODEL_CONE_B,
  COPTMEX_MODEL_CONE_K,
  COPTMEX_MODEL_CONE_Q,

  COPTMEX_MODEL_CONEK_F,
  COPTMEX_MODEL_CONEK_L,
  COPTMEX_MODEL_CONEK_Q,
  COPTMEX_MODEL_CONEK_R,
  COPTMEX_MODEL_CONEK_S,
}; /* model struct fields */


const char* result_fields[] = {
  COPTMEX_RESULT_STATUS,
  COPTMEX_RESULT_SIMITER,
  COPTMEX_RESULT_BARITER,
  COPTMEX_RESULT_NODECNT,
  COPTMEX_RESULT_BESTGAP,
  COPTMEX_RESULT_SOLVETIME,
  COPTMEX_RESULT_OBJVAL,
  COPTMEX_RESULT_BESTBND,
  COPTMEX_RESULT_VARBASIS,
  COPTMEX_RESULT_CONBASIS,
  COPTMEX_RESULT_VALUE,
  COPTMEX_RESULT_REDCOST,
  COPTMEX_RESULT_SLACK,
  COPTMEX_RESULT_DUAL,
  COPTMEX_RESULT_PRIMALRAY,
  COPTMEX_RESULT_DUALFARKAS,

  COPTMEX_RESULT_QCSLACK,

  COPTMEX_RESULT_POOL,
  COPTMEX_RESULT_POOLOBJ,
  COPTMEX_RESULT_POOLXN,

  COPTMEX_RESULT_PSDX,
  COPTMEX_RESULT_PSDRC,
  COPTMEX_RESULT_PSDSLACK,
  COPTMEX_RESULT_PSDPI,
}; /* result struct fields */

const char* iis_fields[] = {
  COPTMEX_IIS_ISMINIIS,
  COPTMEX_IIS_VARLB,
  COPTMEX_IIS_VARUB,
  COPTMEX_IIS_CONSTRLB,
  COPTMEX_IIS_CONSTRUB,
  COPTMEX_IIS_SOS,
  COPTMEX_IIS_INDICATOR,
}; /* IIS result struct fields */

const char* feasibility_fields[] = {
  COPTMEX_FEASRELAX_OBJ,
  COPTMEX_FEASRELAX_LB,
  COPTMEX_FEASRELAX_UB,
  COPTMEX_FEASRELAX_LHS,
  COPTMEX_FEASRELAX_RHS,
}; /* feasibility relaxation result fields */

const char* version_fields[] = {
  COPTMEX_VERSION_MAJOR,
  COPTMEX_VERSION_MINOR,
  COPTMEX_VERSION_TECHNICAL
}; /* version fields */
