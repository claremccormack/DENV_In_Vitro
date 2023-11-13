#ifndef MCMC_H
#define MCMC_H

double Log_Likelihood(double*Obs, double lod, int t, long double mu, double sigma, bool log_scale);
double Log_Prior(const char* distribution, double val, int sample, const char* par1, const char* par2, param_vary& Param_Vary,host& Host);
void Run_MCMC(fit_params Fit, int n_groups, expt*Expt, const vector<string>&v, const vector<string>& d, param_vary& Param_Vary, model_type model, gsl_matrix* Init_Cov);
void Update_Scaling(fit_params Fit, long k, int n, host &Host);
double Check_Boundary(fit_params Fit, const vector<string>& v, int n, double new_val);
void Estimate_Params(fit_params Fit, const vector<string>&v, const vector<string>& d, param_vary& Param_Vary , int n_expt, expt*Expt,  model_type model, gsl_matrix* Init_Cov);
void Update_Block(double size, fit_params Fit, gsl_matrix*Cov, long k, const vector<string>&v,const vector<string>& d, param_vary& Param_Vary , int n_expt, expt* Expt, model_type model);
void Transform_Params(fit_params Fit, long k, double*C, double*N, host &Host, int n, const vector<string> &v);
double Exact_Likelihood(fit_params Fit, long k, int n, host &Host, const vector<string> &v, param_vary& Param_Vary, model_type model);

#endif
