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
