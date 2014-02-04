#include <mex.h>
#include <string.h>
#include <stdio.h>

#define PARS_H "../text1r.h"
#define PARS   text1r_pars_
#define PARS_T text1r_pars_t
#define FIELDS "../text1r.def"

#include PARS_H

/**********************************************************/
/* set a field to matlab/octave structure from an array */
void
set_field(mxArray *mst, int nfld, int N, double *arr){
  int i;
  double *arr1 = mxMalloc(N*sizeof(double));
  mxArray *mxarr = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
  for (i=0; i<N; i++) arr1[i]=arr[i];
  mxSetData(mxarr, arr1);
  mxSetFieldByNumber(mst, 0, nfld, mxarr);
}
/* set a field to matlab/octave structure from an integer array */
void
set_field_i(mxArray *mst, int nfld, int N, int *arr){
  int i;
  double *arr1 = mxMalloc(N*sizeof(double));
  mxArray *mxarr = mxCreateNumericMatrix(N, 1, mxDOUBLE_CLASS, mxREAL);
  for (i=0; i<N; i++) arr1[i]=(double)arr[i];
  mxSetData(mxarr, arr1);
  mxSetFieldByNumber(mst, 0, nfld, mxarr);
}

/* get a field from matlab structure into an integer array */
void
get_field(const mxArray *mst, int nfld, int N, double *arr){
  int i;
  mxArray *v;
  double *arr1;
  v = mxGetFieldByNumber(mst, 0, nfld);
  if (!mxIsDouble(v))
    mexErrMsgTxt("Not a floating point array in a structure field");
  if (mxGetNumberOfElements(v)!=N)
    mexErrMsgTxt("Wrong size of array in a structure field");
  arr1 = mxGetData(v);
  for (i=0;i<N;i++) arr[i] = arr1[i];
}
void
get_field_i(const mxArray *mst, int nfld, int N, int *arr){
  int i;
  mxArray *v;
  double *arr1;
  v = mxGetFieldByNumber(mst, 0, nfld);
  if (!mxIsDouble(v))
    mexErrMsgTxt("Not an array in a structure field");
  if (mxGetNumberOfElements(v)!=N)
    mexErrMsgTxt("Wrong size of array in a structure field");
  arr1 = mxGetData(v);
  for (i=0;i<N;i++) arr[i] = (int)arr1[i];
}

/**********************************************************/
/* Convert c structure to matlab.
   List of field is read from FIELDS file. */
mxArray *
pars_c2mat(PARS_T *cst){
  int i;
  const char *keys[256]; //max number of fields
  i=0;
#define INT(name) keys[i]=#name; i++;
#define DBL(name) keys[i]=#name; i++;
#define ARR(name) keys[i]=#name; i++;
#include FIELDS
#undef INT
#undef DBL
#undef ARR
  mxArray * mst = mxCreateStructMatrix(1,1,i, keys);
  i=0;
#define INT(name)  set_field_i(mst, i, 1, &cst->name); i++;
#define DBL(name)  set_field(mst, i, 1, &cst->name); i++;
#define ARR(name)  set_field(mst, i, cst->n, cst->name); i++;
#include FIELDS
#undef INT
#undef DBL
#undef ARR
  return mst;
}

/* Convert matlab structure to c.
   List of field is read from FIELDS file. */
void
pars_mat2c(PARS_T *cst, const mxArray *mst){
  int i;
  const char * k;

  if (!mxIsStruct(mst))
    mexErrMsgTxt("Structure expected");
  if (mxGetNumberOfElements(mst)!=1)
    mexErrMsgTxt("One element in StructureArray is needed");

  for (i = 0; i < mxGetNumberOfFields(mst); i++){
    k = mxGetFieldNameByNumber(mst, i);
#define INT(name)\
    if (strcmp(k, #name)==0){ get_field_i(mst, i, 1, &cst->name); continue;}
#define DBL(name)\
    if (strcmp(k, #name)==0){ get_field(mst, i, 1, &cst->name); continue;}
#define ARR(name)\
    if (strcmp(k, #name)==0){ get_field(mst, i, cst->n, cst->name); continue;}
#include FIELDS
#undef INT
#undef DBL
#undef ARR
    mexErrMsgTxt("Unknown field");
  }
}

/**********************************************************/
void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[]){

  extern PARS_T PARS;
  int i;

#if FTYPE==1
  // init function
  int j;
  double in[6];
  extern void FUNC(
     double *ttc,double *p,double *nu0,double *r, int *n ,int *itype);

  if (nlhs != 1)
    mexErrMsgTxt("output argument is needed");
  if (nrhs != 6 && nrhs != 5 && nrhs != 4)
    mexErrMsgTxt("Usage: data = text1r_init(ttc, p, nu0, r, n, itype);");

  for (i=0; i<6; i++){
    if (i==4 && nrhs < 5){ in[i]=MAXN; continue; }
    if (i==5 && nrhs < 6){ in[i]=0; continue; }
    if (!mxIsDouble(prhs[i]))
      mexErrMsgTxt("non-numeric argument");
    if (mxGetM(prhs[i])!=1 || mxGetM(prhs[i])!=1)
      mexErrMsgTxt("non-scalar argument");
    in[i] = *mxGetPr(prhs[i]);
  }
  i = (int)in[4]; // double -> int
  j = (int)in[5]; // double -> int
  FUNC(&in[0], &in[1], &in[2], &in[3], &i, &j);
  plhs[0] = pars_c2mat(&PARS);

#elif FTYPE==2

  // set_vortex_* functions
  double in[2];
  extern void FUNC(double *v1, double *v2);
  if (nlhs != 1)
    mexErrMsgTxt("output argument is needed");
  if (nrhs != 3)
    mexErrMsgTxt("wrong number of parameters");
  for (i=1; i<3; i++){
    if (!mxIsDouble(prhs[i]))
      mexErrMsgTxt("non-numeric argument");
    if (mxGetM(prhs[i])!=1 || mxGetM(prhs[i])!=1)
      mexErrMsgTxt("non-scalar argument");
    in[i-1] = *mxGetPr(prhs[i]);
  }
  pars_mat2c(&PARS, prhs[0]);
  FUNC(&in[0], &in[1]);
  plhs[0] = pars_c2mat(&PARS);

#elif FTYPE==3

  // minimize function
  extern void FUNC(int *msglev);
  if (nlhs != 1)
    mexErrMsgTxt("output argument is needed");
  if (nrhs != 1 && nrhs != 2)
    mexErrMsgTxt("wrong number of parameters");
  if (nrhs == 2){
    if (!mxIsDouble(prhs[1]))
      mexErrMsgTxt("non-numeric argument");
    i = *(mxGetPr(prhs[1]));
  }
  else{
    i=-3;
  }

  pars_mat2c(&PARS, prhs[0]);
  FUNC(&i);
  plhs[0] = pars_c2mat(&PARS);

#endif

  return;
}
