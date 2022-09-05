/* SPDX-License-Identifier:  */
#include "coptbex.hpp"


/* workaround */
/**
 * @brief 返回数组中的元素个数.
 * 
 * 如数组维度为 2x5x5，则返回 30
 */
size_t bxGetNumberOfElements(const bxArray *ba) {
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
    constexpr size_t TOTAL_PLUGIN_FUNCTIONS = 1;
    bexfun_info_t* func_list_dyn = new bexfun_info_t[TOTAL_PLUGIN_FUNCTIONS + 1];

    size_t i = 0;
    func_list_dyn[i].name = "copt_version";
    func_list_dyn[i].ptr  = copt_version;
    func_list_dyn[i].help = copt_version_help;

    // 最后一个元素, `name` 字段必须为空字符串 `""`
    i++;
    func_list_dyn[i].name = "";
    func_list_dyn[i].ptr  = nullptr;
    func_list_dyn[i].help = nullptr;

    assert((TOTAL_PLUGIN_FUNCTIONS == i));
    return func_list_dyn;
} /* bxPluginFunctions */
