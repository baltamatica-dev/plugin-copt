/* SPDX-License-Identifier:  */
#pragma once
#include <iostream>
#include <cassert>
#include <cmath>
#include <bex/bex.hpp>


#define BUILD_WITH_BEX_WARPPER
#define bxLogical   bool

#define bxFree      free
#define bxCalloc    calloc
#define bxIsNaN     std::isnan


/* workaround */

extern void mexErrMsgIdAndTxt(const char *errorid, const char *errormsg, ...);
extern const char* bxGetFieldNameByNumber(const bxArray *pm, int fieldnumber);
extern int bxGetNumberOfElements(const bxArray *pm);
extern int bxAddField(bxArray *pm, const char *fieldname);
extern bool bxIsScalar(const bxArray*);
extern double bxGetScalar(const bxArray*);
extern int mexEvalString(const char *command);


/* 导出函数 */

BALTAM_PLUGIN_FCN(copt_computeiis);
BALTAM_PLUGIN_FCN(copt_feasrelax);
BALTAM_PLUGIN_FCN(copt_read);
BALTAM_PLUGIN_FCN(copt_solve);
BALTAM_PLUGIN_FCN(copt_write);


#define MODEL_FIELD_NUM  48
extern const char* model_fields[MODEL_FIELD_NUM];
#define RESULT_FIELD_NUM  24
extern const char* result_fields[RESULT_FIELD_NUM];
#define IIS_FIELD_NUM  7
extern const char* iis_fields[IIS_FIELD_NUM];
#define FEASRELAX_FIELD_NUM  5
extern const char* feasibility_fields[FEASRELAX_FIELD_NUM];
#define VERSION_FIELD_NUM  3
extern const char* version_fields[VERSION_FIELD_NUM];
