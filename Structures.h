#ifndef STRUCTURES_H
#define STRUCTURES_H

#define RUNS 1
#define RUN_IN 1
#define DAYS 10+RUN_IN
#define DT .04166667 /*timestep*/
#define STEPS 24*(DAYS)+2

#include <vector>
#include <string>

using namespace std;

typedef struct MODEL_TYPE {

	int in_vitro_model; /*0-with latent period, 1-no latent period*/
	int drug_moa; /*allow for different mechanisms of drug action*/

	int n_expt; /*number of experiments*/
	
}model_type;

typedef struct MODEL_PARAMS {

		double s_n; /*cell growth rate*/
		double kappa_t; /*target cell mortality rate*/
		double K; /*carrying capacity of cell dish*/
		
		double beta; /*target cell infection rate (per virion)*/
		double kappa_i; /*increase in mortality due to infection*/
		double tau; /*inverse of latent period*/

		double omega; /*intracellular virus proliferation rate*/
		double epsilon; /*duration of time spent in intracellular phase*/
		
		double prop_v; /*proportion of intracellular virus which becomes extracellular*/
		double kappa_v; /*extracellular virus clearance rate*/		
		
		double R0;

}model_params;

typedef struct PARAM_VARY {

	std::vector<std::string>Name;
	std::vector<std::string>Dist;
	int*logscale;
	int no_params;

}param_vary;

typedef struct DOSE {

	int dose_no;
	int dose_time;
	double dose_size;

}dose;

typedef struct PK {
	
	double***A1;

	/*parameters for concentration-effect model*/

	double hill_coeff;
	double ec50;
	double emax;

}pk;

typedef struct HOST_DATA {

	int no_measurement_times;
	int Time[STEPS];

	std::vector<std::string>Strain;
	std::vector<std::string>Refresh;

	/*limits of detection and measurement error levels*/

		double lod_v; /*limit of detection (extracellular virus)*/
		double lod_x; /*limit of detection (intracellular virus)*/
	
		double error_v; /*measurement error (extracellular virus)*/
		double error_x; /*measurement error (intracellular virus)*/
		
	/*arrays to store data*/

		double X[STEPS]; /*intracellular virus*/
		double V[STEPS]; /*extracellular virus*/
		double D[STEPS]; /*drug concentration*/

}host_data;

typedef struct HOST_AVG {

	/*arrays to store average values*/

	long double** T; /*target cells*/
	long double** E; /*exposed cells*/
	long double** I; /*infectious cells*/
	long double** X; /*intracellular virus*/
	long double** V; /*extracellular virus*/
	long double** D; /*drug*/


}host_avg;

typedef struct HOST_PREDICT {

	/*stores predicted values*/

	long double*T; /*target cells*/
	long double*E; /*exposed cells*/
	long double*I; /*infectious cells*/
	long double*X; /*intracellular virus (in-vitro model)*/
	long double*V; /*virus*/
	long double*D; /*drug*/
	long double*R;
	double*Eff;

}host_predict;

typedef struct PARAM_CHAIN {

	double*curr_val; /*current parameter value*/
	double*new_val; /*proposed parameter value*/
	double*curr_ll; /*log-likelihood value based of current parameter values*/
	double*new_ll; /*log-likelihood value based of proposed parameter values*/
	double*curr_pr; /*prior value based of current parameter values*/
	double*new_pr; /*prior value based of proposed parameter values*/
	double*curr_posterior; /*posterior value based of current parameter values*/
	double*new_posterior; /*posterior value based on proposed parameter values*/
	double*curr_sd; /*current jump size*/
	double*new_sd; /*new jump size*/
	double*no_acpt; /*number of accepted values*/
	double curr_R0; /*R0 value based on current parameter values*/
	double new_R0; /*R0 value based on proposed parameter values*/

	double**C; /*chain of estimated parameter values*/
	double*LC; /*chain of log-likelihood values*/
	double*R0; /*chain of R0 values*/
	double**Sub_C; /*subset of estimated parameter values*/
	double*Sub_LC; /*subset of log-likelihood value*/
	double*Sub_R0; /*subset of R0 values*/
	double**SD; /*chain of jump sizes*/
	double**Posterior; /*posterior values*/
	double**Prop_Acpt; /*proportion of proposed parameter values accepted*/

	/*credible intervals*/

		double*Lower; /*lower bound of 95% credible interval for estimated parameter values*/
		double*Median; /*median estimated parameter values*/
		double*Upper; /*upper bound of 95% credible interval for estimated parameter values*/
		double LL_Median; /*median estimated log-likelihood value*/
		double R0_Median; /*median estimated R0 value*/
		double LL_Lower; /*lower bound of 95% credible interval for log-likelihood values*/
		double R0_Lower; /*lower bound of 95% credible interval for R0 values*/
		double LL_Upper; /*upper bound of 95% credible interval for log-likelihood values*/
		double R0_Upper; /*upper bound of 95% credible interval for R0 values*/

	double** W;
	double** Param_Sample;
	double* Mean_Param_Sample;

}param_chain;

typedef struct FIT_PARAMS {

	int no_fit; /*number of parameters we want to estimate*/
	double iterations; /*number of MCMC iterations*/
	double burn_in; /*length of burn-in for each MCMC chain*/
	int no_post_samples; /*number of samples from posterior parameter distributions*/
	double req_acpt_rate; /*desired MCMC accpetance rate*/
	int block_update; /*0-sequentially update parameter values, 1-block update parameters*/
	int thin; /*factor by which we thin MCMC chains*/

	double*start_values; /*initial value of each parameter we're estimating*/
	int*logscale; /*0- fit parameter on linear scale, 1- fit parameter on log10 scale*/
	int*level; /*0-fitting parameter globally, 1- fitting parameter locally to each experiment, 2-fitting parameter locally to each group, within an experiment, 3- fitting parameter to each individual host;*/
	int* block;
	vector<vector<string>>lower; /*lower bound for prior distribution*/
	vector<vector<string>>upper; /*upper bound for prior distribution*/
	int*sample;
	double*init_sd; /*initial jump size for each parameter we're estimating*/
	double*max_sd; /*max jump size for each parameter we're estimating*/
	double*min_sd; /*min jump size for each parameter we're estimating*/
	double scaling; /*controls how quickly we adjust jump size - large values give smaller, slower adjustments*/

	int no_particles;
	int trajectory;
	int cov_start;
	int cov_stop;
	int cov_update;

}fit_params;

/*hierarchical structure: experiment-group-individual host*/

typedef struct HOST {

	std::string strain;
	std::string refresh;

	model_params Params; /*viral dynamics model parameters*/
	pk PK; /*pk model parameters*/
	host_data Data; /*data*/

	param_chain Par_Chain; /*results of MCMC*/
	host_predict*Simulation; /*example simulations using parameter values estimated during MCMC*/
	host_predict*Bounds; /*stores posterior predictive intervals*/

	/*basic info about each host*/
	
		double viral_inoculum;;
		int time_infection;

		int no_doses;
		int first_dose;
		int last_dose;
		double dose_size;
		double dose_frequency;
		double dosage[STEPS];	
		dose* Doses;
		
	/*arrays to store cell populations for each host*/

		long double***T; /*target cells*/
		long double***E; /*exposed cells*/
		long double***I; /*infectious cells*/
		long double***X; /*intracellular virus*/
		long double***V; /*extracellular virus*/
		long double***D; /*drug*/
	
	double** LL;

}host;

typedef struct POP {

	int no_hosts; /*number of hosts within a group (within an experiment)*/
	host* Hosts;
	host_avg* Hosts_Avg;
	host_avg Avg;

}pop;

typedef struct EXPT {

	int no_groups; /*number of groups within an experiment*/
	pop* Pop;
	
}expt;

#endif
