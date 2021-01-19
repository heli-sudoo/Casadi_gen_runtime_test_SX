/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

int Dyn_BS(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int Dyn_BS_alloc_mem(void);
int Dyn_BS_init_mem(int mem);
void Dyn_BS_free_mem(int mem);
int Dyn_BS_checkout(void);
void Dyn_BS_release(int mem);
void Dyn_BS_incref(void);
void Dyn_BS_decref(void);
casadi_int Dyn_BS_n_out(void);
casadi_int Dyn_BS_n_in(void);
casadi_real Dyn_BS_default_in(casadi_int i);
const char* Dyn_BS_name_in(casadi_int i);
const char* Dyn_BS_name_out(casadi_int i);
const casadi_int* Dyn_BS_sparsity_in(casadi_int i);
const casadi_int* Dyn_BS_sparsity_out(casadi_int i);
int Dyn_BS_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif