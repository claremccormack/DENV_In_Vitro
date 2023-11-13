#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cmath>
#include<random>
#include<algorithm>
#include<time.h>
#include<fstream>
#include<string>
#include<sstream>
#include<direct.h>
#include<windows.h>
#include <omp.h>
#include "randlib_par.h"
#include "Structures.h"
#include "Functions.h"
#include "In_Vitro_Model.h"
#include "MCMC.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_math.h"

random_device rd2;
mt19937 gen2(rd2());

gsl_ran_discrete_t * F;
gsl_rng * g1 = gsl_rng_alloc(gsl_rng_mt19937);

uniform_real_distribution<double> std_uniform(0, 1);

double pi = 4 * atan(1);

void Run_MCMC(fit_params Fit, int n_expt, expt*Expt, const vector<string>&v, const vector<string>& d, param_vary& Param_Vary, model_type model, gsl_matrix* Init_Cov) {

	gsl_rng_set(g1, long(time(NULL)));

	for (int i = 0; i < n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				for (int n = 0; n < Fit.no_fit; n++) {
					
					Expt[i].Pop[j].Hosts[k].Par_Chain.no_acpt[n] = 0;
				}
			}
		}
	}

	Estimate_Params(Fit, v, d, Param_Vary, n_expt, Expt, model, Init_Cov);
}

void Estimate_Params(fit_params Fit, const vector<string>& v, const vector<string>& d, param_vary& Param_Vary, int n_expt, expt* Expt, model_type model, gsl_matrix* Init_Cov) {

	stringstream ss_n, ss_m;
	int k0;
	int k_step = int(Fit.iterations / Fit.thin);

	/*calculate how many parameters in the block*/

		double size = Block_Size(Expt, Fit, model);

	/*allocate memory for covariance matrix*/

			gsl_matrix* Cov,*I,*Corr;
			gsl_vector* mean;

			Cov = gsl_matrix_alloc(size, size);/*covariance matrix*/
			Corr = gsl_matrix_alloc(size, size);/*correlation matrix*/
			I = gsl_matrix_alloc(size, size);/*identity matrix*/

			mean = gsl_vector_alloc(size);

			gsl_matrix_set_identity(I);/*set identity matrix*/

		/*Until iteration 'cov_start', update parameter values using identity matrix*/

			for (long k = 0; k < Fit.cov_start+1; k++) {

				cout << k << endl;

				Update_Block(size, Fit, Init_Cov, k, v, d, Param_Vary, n_expt, Expt, model);

			}

			/*calculate covariance matrix*/

				Covariance_Matrix(size, 0, Fit.cov_start/double(Fit.thin),Expt,Fit,model,Cov,Corr);


			/*From iteration 'cov_start' until iteration 'cov_stop' update parameters using covariance matrix, which is been updated every 'cov_update' iterations*/

			int k_max;

			if (Fit.cov_start > Fit.iterations) {k_max = int(Fit.iterations - 1);}

			else {k_max = int(Fit.cov_stop);}

			for (int k = Fit.cov_start+1; k < k_max; k++) {

				cout << k << endl;

				Update_Block(size,Fit, Corr, k, v, d, Param_Vary, n_expt, Expt, model);

				/*Every 'cov_update' iterations, calculate covariance matrix using a maximum of the previous 200,000 iterations*/

					int k0;

					if (k != (Fit.iterations - 1) && (k % Fit.cov_update == 0)) {

						if (k <= 200000) { k0 = 0; }

						else { k0 = (k - 200000)/double(Fit.thin); }

						Covariance_Matrix(size, k0, k/double(Fit.thin), Expt, Fit, model, Cov, Corr);

					}
			}

			/*For remaining iterations, fix covariance matrix at most recent value*/

			for (int k = Fit.cov_stop; k < Fit.iterations; k++) {

				cout << k << endl;

				Update_Block(size, Fit, Corr, k, v, d, Param_Vary, n_expt, Expt, model);

			}

			ss_n << n_expt << ',' << v[0] << ',' << ',' << Fit.level[0] << ',' << Fit.no_fit << ',' << Fit.start_values[0] << ',' << Fit.init_sd[0];

			/*delete memory stored*/

			gsl_matrix_free(Cov);
			gsl_matrix_free(Corr);
			gsl_matrix_free(I);
			gsl_vector_free(mean);
		}

void Update_Block(double size, fit_params Fit, gsl_matrix*Cov, long k, const vector<string>&v, const vector<string>& d, param_vary& Param_Vary, int n_expt, expt* Expt, model_type model) {

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///For method see pages 9-11 of "Adaptive Optimal Scaling of Metropolis-Hastings Algorithms Using the Robbins-Monro Process (Garthwaite et. al, 2010)*///
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double epsilon, ratio_post, correction, p_accept, temp, new_val, x, curr_val,sum_ll, sum_pr;
	int k0, thread_no;

	//int thread_count = omp_get_max_threads();
	int thread_count=1;
	omp_set_num_threads(thread_count);

	gsl_vector *mean, *jump;
	gsl_matrix *I, *Cov_S;

	Cov_S = gsl_matrix_alloc(size, size);/*stores scaled covariance matrix*/
	gsl_matrix_memcpy(Cov_S, Cov); /*copy entries from Cov to Cov_S*/

	I = gsl_matrix_alloc(size, size);
	gsl_matrix_set_identity(I);/*identity matrix*/

	jump = gsl_vector_alloc(size); /*stores jump size for each parameter*/

	mean = gsl_vector_alloc(size); /*stores mean of proposal distribution*/

	for (int s = 0; s < size; s++) { gsl_vector_set(mean, s, 0); }/*set mean of proposal distribution to zero*/

	/*tune jump size*/

	Update_Scaling(Fit, k, 0, Expt[0].Pop[0].Hosts[0]);

	/*add a small amount (epsilon) to diagonal elements of the scaled covariance matrix to ensure it's always positive definite (otherwise can't take Cholesky decomposition of the matrix, needed in next step)*/

	epsilon = pow(0.001, 2) / (k + 1);

	gsl_matrix_scale(I, epsilon);

	gsl_matrix_add(Cov_S, I);

/*scale covariance matrix*/

	gsl_matrix_scale(Cov_S, pow(Expt[0].Pop[0].Hosts[0].Par_Chain.new_sd[0], 2));

	/*get Cholesky decomposition of the covariance matrix and draw jump sizes from multivariate normal distribution with covariance matrix Cov_S with mean 0 */

		gsl_linalg_cholesky_decomp(Cov_S);

		gsl_ran_multivariate_gaussian(g1, mean, Cov_S, jump);

	/*proposed new parameter values*/

		int s = 0;

		for (int n = 0; n < Fit.no_fit; n++) {

			switch (Fit.level[n]) {
			 
			case 0: {

				curr_val = Expt[0].Pop[0].Hosts[0].Par_Chain.curr_val[n];
				x = gsl_vector_get(jump, s); /*get jump size corresponding to index n*/
				new_val = curr_val + x;

				new_val = Check_Boundary(Fit, v, n, new_val);

				for (int i = 0; i < n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int l = 0; l < Expt[i].Pop[j].no_hosts;l++) {

							Expt[i].Pop[j].Hosts[l].Par_Chain.new_val[n] = new_val;
							Expt[i].Pop[j].Hosts[l].Par_Chain.curr_sd[n] = Expt[0].Pop[0].Hosts[0].Par_Chain.new_sd[0];
						}
					}
				}

				s = s + 1;

				break;
			}

			case 1: {

				for (int i = 0; i < n_expt;i++) {

				curr_val = Expt[i].Pop[0].Hosts[0].Par_Chain.curr_val[n];
				x = gsl_vector_get(jump, s); /*get jump size corresponding to index n*/
				new_val = curr_val + x;

				new_val = Check_Boundary(Fit, v, n, new_val);
	
		
					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int l = 0; l < Expt[i].Pop[j].no_hosts;l++) {

							Expt[i].Pop[j].Hosts[l].Par_Chain.new_val[n] = new_val;
							Expt[i].Pop[j].Hosts[l].Par_Chain.curr_sd[n] = Expt[0].Pop[0].Hosts[0].Par_Chain.new_sd[0];
						}
					}

					s = s + 1;
				}

				break;
			}

			case 2: {

				for (int i = 0; i < n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						curr_val = Expt[i].Pop[j].Hosts[0].Par_Chain.curr_val[n];
						x = gsl_vector_get(jump, s); /*get jump size corresponding to index n*/
						new_val = curr_val + x;
						
						new_val = Check_Boundary(Fit, v, n, new_val);

						for (int l = 0; l < Expt[i].Pop[j].no_hosts;l++) {

							Expt[i].Pop[j].Hosts[l].Par_Chain.new_val[n] = new_val;
							Expt[i].Pop[j].Hosts[l].Par_Chain.curr_sd[n] = Expt[0].Pop[0].Hosts[0].Par_Chain.new_sd[0];
						}

						s = s + 1;
					}
				}

				break;
			}
	
			}	
		}

	/*run MCMC algorithm*/
	
		if (k == 0) {
			
				for (int i = 0; i < n_expt;i += 1) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					#pragma omp parallel for private(thread_no) schedule(static,1)

					for (int thread_no = 0; thread_no < thread_count; thread_no++) {

					for (int l = thread_no; l < Expt[i].Pop[j].no_hosts; l += thread_count) {

							Transform_Params(Fit, k, Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val, Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val, Expt[i].Pop[j].Hosts[l], 0, v);

							Expt[i].Pop[j].Hosts[l].Par_Chain.curr_ll[0]=Exact_Likelihood(Fit, k, 0, Expt[i].Pop[j].Hosts[l], v, Param_Vary, model);

							Expt[i].Pop[j].Hosts[l].Par_Chain.curr_R0 = Calculate_Rt(Expt[i].Pop[j].Hosts[l], Expt[i].Pop[j].Hosts[l].time_infection+24*RUN_IN, model);

							for (int m = 0; m < Fit.no_fit;m++) {

								Expt[i].Pop[j].Hosts[l].Par_Chain.curr_pr[m] = 0;

								if (Fit.sample[m] == 0) {
									
									Expt[i].Pop[j].Hosts[l].Par_Chain.curr_pr[m]+=Log_Prior(d[m].c_str(), Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val[m], Fit.sample[m], Fit.lower[m][i].c_str(), Fit.upper[m][i].c_str(), Param_Vary, Expt[i].Pop[j].Hosts[l]); }
								 
							}

						}
					}
				}
			}

			sum_ll = 0;
			sum_pr = 0;

			for (int i = 0; i < n_expt;i++) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

						sum_ll += Expt[i].Pop[j].Hosts[l].Par_Chain.curr_ll[0];

						for (int n = 0; n < Fit.no_fit;n++) {

							sum_pr += Expt[i].Pop[j].Hosts[l].Par_Chain.curr_pr[n];
						}
					}
				}
			}

		for (int i = 0; i < n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

					for (int n = 0; n < Fit.no_fit;n++) {

						Expt[i].Pop[j].Hosts[l].Par_Chain.curr_ll[n] = sum_ll;
						Expt[i].Pop[j].Hosts[l].Par_Chain.LC[0] = sum_ll;
						Expt[i].Pop[j].Hosts[l].Par_Chain.R0[0] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_R0;
						Expt[i].Pop[j].Hosts[l].Par_Chain.curr_posterior[n] = sum_ll + sum_pr;
		
					}
				}
			}
		}

	}

		else {

			for (int i = 0; i < n_expt; i += 1) {

				for (int j = 0; j < Expt[i].no_groups; j++) {

					#pragma omp parallel for private(thread_no) schedule(static,1)

					for (int thread_no = 0; thread_no < thread_count; thread_no++) {

						for (int l = thread_no; l < Expt[i].Pop[j].no_hosts; l += thread_count) {

						Transform_Params(Fit, k, Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val, Expt[i].Pop[j].Hosts[l].Par_Chain.new_val, Expt[i].Pop[j].Hosts[l], 0, v);

						Expt[i].Pop[j].Hosts[l].Par_Chain.new_ll[0] = Exact_Likelihood(Fit, k + 1, 0, Expt[i].Pop[j].Hosts[l], v, Param_Vary, model);

						Expt[i].Pop[j].Hosts[l].Par_Chain.new_R0 = Calculate_Rt(Expt[i].Pop[j].Hosts[l], Expt[i].Pop[j].Hosts[l].time_infection + 24 * RUN_IN, model);

							for (int m = 0; m < Fit.no_fit; m++) {

								Expt[i].Pop[j].Hosts[l].Par_Chain.new_pr[m] = 0;

								if (Fit.sample[m] == 0) {

									Expt[i].Pop[j].Hosts[l].Par_Chain.new_pr[m] += Log_Prior(d[m].c_str(), Expt[i].Pop[j].Hosts[l].Par_Chain.new_val[m], Fit.sample[m], Fit.lower[m][i].c_str(), Fit.upper[m][i].c_str(), Param_Vary, Expt[i].Pop[j].Hosts[l]);

								}
							}

						}
					}
				}
			}

			sum_ll = 0;
			sum_pr = 0;

			for (int i = 0; i < n_expt; i++) {

				for (int j = 0; j < Expt[i].no_groups; j++) {

					for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

						sum_ll += Expt[i].Pop[j].Hosts[l].Par_Chain.new_ll[0];

						for (int n = 0; n < Fit.no_fit; n++) {

							sum_pr += Expt[i].Pop[j].Hosts[l].Par_Chain.new_pr[n];
						}
					}
				}
			}

			for (int i = 0; i < n_expt; i++) {

				for (int j = 0; j < Expt[i].no_groups; j++) {

					for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

						for (int n = 0; n < Fit.no_fit; n++) {

							Expt[i].Pop[j].Hosts[l].Par_Chain.new_ll[n] = sum_ll;
							Expt[i].Pop[j].Hosts[l].Par_Chain.new_posterior[n] = sum_ll + sum_pr;
						}
					}
				}

			}

			/*calculate probability of accepting new value*/

			ratio_post = Expt[0].Pop[0].Hosts[0].Par_Chain.new_posterior[0] - Expt[0].Pop[0].Hosts[0].Par_Chain.curr_posterior[0];
			correction = 0;
			p_accept = ratio_post + correction;

			if (p_accept > 0) { p_accept = 0; }

			temp = log(std_uniform(gen2));/*working on log scale*/

			if (temp < p_accept) {

				for (int i = 0; i < n_expt; i++) {

					for (int j = 0; j < Expt[i].no_groups; j++) {

						for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

							for (int n = 0; n < Fit.no_fit; n++) {

								Expt[i].Pop[j].Hosts[l].Par_Chain.no_acpt[n] = Expt[i].Pop[j].Hosts[l].Par_Chain.no_acpt[n] + 1; /*increase no. of acceptances*/
								Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val[n] = Expt[i].Pop[j].Hosts[l].Par_Chain.new_val[n];
								Expt[i].Pop[j].Hosts[l].Par_Chain.curr_posterior[n] = Expt[i].Pop[j].Hosts[l].Par_Chain.new_posterior[n];
								Expt[i].Pop[j].Hosts[l].Par_Chain.curr_ll[n] = Expt[i].Pop[j].Hosts[l].Par_Chain.new_ll[n];
								Expt[i].Pop[j].Hosts[l].Par_Chain.curr_R0 = Expt[i].Pop[j].Hosts[l].Par_Chain.new_R0;

								if (k % Fit.thin == 0) {

									k0 = int(k / Fit.thin);

									Expt[i].Pop[j].Hosts[l].Par_Chain.LC[k0] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_ll[n];
									Expt[i].Pop[j].Hosts[l].Par_Chain.R0[k0] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_R0;
									Expt[i].Pop[j].Hosts[l].Par_Chain.C[k0][n] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val[n];
									Expt[i].Pop[j].Hosts[l].Par_Chain.SD[k0][n] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_sd[n];
									Expt[i].Pop[j].Hosts[l].Par_Chain.Prop_Acpt[k0][n] = (Expt[i].Pop[j].Hosts[l].Par_Chain.no_acpt[n]) / (k + 1);

								}

							}
						}
					}
				}
			}

			else {

				for (int i = 0; i < n_expt; i++) {

					for (int j = 0; j < Expt[i].no_groups; j++) {

						for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

							for (int n = 0; n < Fit.no_fit; n++) {

								if (k % Fit.thin == 0) {

									k0 = int(k / Fit.thin);

									Expt[i].Pop[j].Hosts[l].Par_Chain.LC[k0] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_ll[n];
									Expt[i].Pop[j].Hosts[l].Par_Chain.R0[k0] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_R0;
									Expt[i].Pop[j].Hosts[l].Par_Chain.C[k0][n] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_val[n];
									Expt[i].Pop[j].Hosts[l].Par_Chain.SD[k0][n] = Expt[i].Pop[j].Hosts[l].Par_Chain.curr_sd[n];
									Expt[i].Pop[j].Hosts[l].Par_Chain.Prop_Acpt[k0][n] = (Expt[i].Pop[j].Hosts[l].Par_Chain.no_acpt[n]) / (k + 1);

								}
							}
						}
					}
				}
			}
		}

	/*delete memory*/

	gsl_vector_free(mean);
	gsl_vector_free(jump);
	gsl_matrix_free(I);
	gsl_matrix_free(Cov_S);
}

double Log_Likelihood(double* Obs, double LOD, int i, long double mu, double sigma, bool log_scale) {

	double lik, log_lik, l;

	if (Obs[i] == 0) { Obs[i] = 1e-6; }

	if (mu == 0) { mu = 1e-6; }

	if (log_scale == 1) {

		if (Obs[i] <= log10(LOD)) {

			l = 0.5 * erfc(-(log10(LOD) - log10(mu)) / (sigma * sqrt(2)));

			lik = l;

			if ((1 - l) < 0.00001) { lik = 1; }

			if (lik > 0) { log_lik = log(lik); }

			else { log_lik = -200000; }

		}

		else {

			log_lik = -0.5 * log(2 * pi) - 0.5 * log(pow(sigma, 2)) - 0.5 * (pow(((Obs[i] - log10(mu)) / sigma), 2));
		}
	}

	else {

		if (Obs[i] <= LOD) {

			l = 0.5 * erfc(-(LOD - mu) / (sigma * sqrt(2)));

			lik = l;

			if ((1 - l) < 0.00001) { lik = 1; }

			if (lik > 0) { log_lik = log(lik); }

			else { log_lik = -200000; }

		}

		else {

			log_lik = -0.5 * log(2 * pi) - 0.5 * log(pow(sigma, 2)) - 0.5 * (pow(((Obs[i] - mu) / sigma), 2));
		}
	}

	return log_lik;

}

double Check_Boundary(fit_params Fit, const vector<string>& v, int n, double new_val) {

	vector<string> Params;
	double temp;

	Params = {"prop_v"}; /*parameters in [0,1] range*/

	if (find(Params.begin(), Params.end(), v[n]) != Params.end()) {

		if (Fit.logscale[n] == 0) {

			if (new_val < 0) { temp = 1 + new_val; }

			else if (new_val > 1) { temp = new_val - 1; }

			else { temp = new_val; }

			while (temp < 0 || temp >1) {

				if (temp < 0) { temp = 1 + temp; }

				else { temp = temp - 1; }
			}
		}

		else {

			if (new_val < -25) { temp = new_val + 25; }

			else if (new_val > 0) { temp = -25 + new_val; }

			else { temp = new_val; }

			while (temp < -25 || temp >0) {

				if (temp < -25) { temp = temp - 25; }

				else { temp = -25 + temp; }

			}
		}
	}

	else {

		/*ensure parmeter value is positive*/

		if (Fit.logscale[n] == 0) {

			if (new_val < 0) { temp = -new_val; }

			else temp = new_val;
		}

		else {

			if (new_val < -25) { temp = -new_val; }

			else temp = new_val;
		}

		temp = new_val;

	}

	return temp;
}

void Transform_Params(fit_params Fit, long k, double*C, double*N, host &Host, int n, const vector<string> &v) {

		for (int i = 0; i < Fit.no_fit; i++) {

			if (strcmp(v[i].c_str(), "beta") == 0) { if (Fit.logscale[i] == 0) { Host.Params.beta = N[i]; } else { Host.Params.beta = pow(10, N[i]);}}
			if (strcmp(v[i].c_str(), "hill_coeff") == 0) { if (Fit.logscale[i] == 0) { Host.PK.hill_coeff = N[i]; } else { Host.PK.hill_coeff = pow(10, N[i]); } }
			if (strcmp(v[i].c_str(), "ec50") == 0) { if (Fit.logscale[i] == 0) { Host.PK.ec50 = N[i]; } else { Host.PK.ec50= pow(10, N[i]); } }
			if (strcmp(v[i].c_str(), "omega") == 0) { if (Fit.logscale[i] == 0) { Host.Params.omega = N[i]; } else { Host.Params.omega = pow(10, N[i]); } }
			if (strcmp(v[i].c_str(), "prop_v") == 0) { if (Fit.logscale[i] == 0) { Host.Params.prop_v = N[i]; } else { Host.Params.prop_v = pow(10, N[i]); } }
			if (strcmp(v[i].c_str(), "error_v") == 0) { if (Fit.logscale[i] == 0) { Host.Data.error_v = N[i]; } else { Host.Data.error_v = pow(10, N[i]); } }
			if (strcmp(v[i].c_str(), "error_x") == 0) { if (Fit.logscale[i] == 0) { Host.Data.error_x = N[i]; } else { Host.Data.error_x = pow(10, N[i]); } }
		
		}

}

void Update_Scaling(fit_params Fit, long k, int n, host &Host) {

	/*Update jump sizes using Robbins-Monro process algorithms in "Adaptive Optimal Scaling of Metropolis-Hastings Algorithms Using the Robbins-Monro Process (Garthwaite et. al, 2010)*/

	double c, q1, q2, q3, T;

	T = Fit.scaling; /*larger values of T reduce the size of each adjustment*/

	if (k < 2) { Host.Par_Chain.new_sd[n] = Host.Par_Chain.curr_sd[n]; }

		else {

			q1 = -gsl_cdf_ugaussian_Pinv(Fit.req_acpt_rate / 2);

			q2 = (sqrt(2 * 22 / 7)*exp(pow(q1, 2) / 2)) / (2 * q1);

			q3 = 1 / double(Fit.no_fit*Fit.req_acpt_rate*(1 - Fit.req_acpt_rate));

			c = Host.Par_Chain.curr_sd[n] * ((1 - (1 / double(Fit.no_fit)))*q2 + q3);

			//cout << q1 << ',' << q2 << ',' << q3 << ',' << c << endl;

			/*If previous proposed value was accepted, increase variance*/
			if (Host.Par_Chain.curr_val[n] == Host.Par_Chain.new_val[n]) { Host.Par_Chain.new_sd[n] = Host.Par_Chain.curr_sd[n] + (c*(1 - Fit.req_acpt_rate) / max(T, (k - 1) / double(Fit.no_fit))); }

			/*If previous proposed value was rejected, decrease variance*/
			else { Host.Par_Chain.new_sd[n] = Host.Par_Chain.curr_sd[n] - (c*(Fit.req_acpt_rate) / max(T, (k - 1) / double(Fit.no_fit))); }

			/*Do not let jump size go beyond specified min/max values*/

			if (Host.Par_Chain.new_sd[n] < Fit.min_sd[n]) { Host.Par_Chain.new_sd[n] = Fit.min_sd[n]; }

			if (Host.Par_Chain.new_sd[n] > Fit.max_sd[n]) { Host.Par_Chain.new_sd[n] = Fit.max_sd[n]; }
	}

}

double Exact_Likelihood(fit_params Fit, long k, int n, host &Host, const vector<string> &v, param_vary& Param_Vary, model_type model) {

	/*get initial values for each compartment based on proposed parameter value*/

		Init_Values(Host, Fit.no_particles,model);

		Drug_Dosage(Host,model);

	/*run model*/

		In_Vitro_Model(Host, 0, 0, STEPS - 2, 1, model);

	/*calculate log-likelihood*/

		double log_lik = 0;

		for (int i = 0; i < Host.Data.no_measurement_times; i++) {

			int t = Host.Data.Time[i] + 24*RUN_IN;
		
			if (isnan(Host.Data.V[i]) == 0) {

				log_lik+=Log_Likelihood(Host.Data.V, Host.Data.lod_v, i, Host.V[0][t][0], Host.Data.error_v,1);
			}
		
			if (isnan(Host.Data.X[i]) == 0) {

				log_lik += Log_Likelihood(Host.Data.X, Host.Data.lod_x, i, Host.X[0][t][0], Host.Data.error_x,1);
			}

		}

	return log_lik;

}

double Log_Prior(const char*distribution, double val, int sample, const char* par1, const char* par2, param_vary& Param_Vary, host& Host) {

	double result,a,b;

	result = 0;

		a=stod(par1);
		b=stod(par2);

		if (strcmp(distribution, "Uniform") == 0) { result = gsl_ran_flat_pdf(val, a, b); }
		if (strcmp(distribution, "Normal") == 0) { result = gsl_ran_gaussian_pdf(val - a, b); }
		
	if (result == 0) { result = -50000; }

	else { result = log(result); }

	return result;

}