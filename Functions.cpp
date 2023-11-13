#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<cmath>
#include<ctime>
#include<random>
#include<algorithm>
#include<string>
#include<fstream>
#include"randlib_par.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_math.h"
#include <gsl/gsl_statistics.h>
#include"Structures.h"
#include"Functions.h"
#include"In_Vitro_Model.h"
#include"MCMC.h"

random_device rd1;
mt19937 gen1(rd1());

gsl_ran_discrete_t * F1;
gsl_rng * g2 = gsl_rng_alloc(gsl_rng_mt19937);

const double na = sqrt(-1);

void Read_Data(string filename, expt*& Expt, int& n_expt) {

	ifstream Data_File;

	Data_File.open(filename);

	string line, experiment, group, individual, dose, strain, refresh,inoculum, time_infection, time_measurement, v, x, d, lod_v,lod_x;

	int i, j, k, prev_i, prev_j, prev_k;

	int t = 0;

	prev_i = 0;
	prev_j = 0;
	prev_k = 0;

	n_expt = 0;

	/*count how many experiments*/

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);

			if (i == 0) { n_expt = 1; }

			if (i > prev_i) {

				n_expt += 1;
			}

			prev_i = i;
		}

		Data_File.close();

		Expt = new expt[n_expt];

	/*count how many groups within each experiment*/

		prev_i = 0;

		Data_File.open(filename);

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);
			j = stoi(group);

			if (j == 0) { Expt[i].no_groups = 1; }

			if (i == prev_i && j > prev_j) {

				Expt[i].no_groups += 1;
			}

			prev_i = i;
			prev_j = j;
		}

		Data_File.close();

		for (int i = 0; i < n_expt; i++) { Expt[i].Pop = new pop[Expt[i].no_groups];}


	/*count how many individuals within each group*/

		prev_i = 0;
		prev_j = 0;

		Data_File.open(filename);

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\t');
			getline(Data_File, individual, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);
			j = stoi(group);
			k = stoi(individual);

			if (k == 0) { Expt[i].Pop[j].no_hosts = 1; }

			if (j == prev_j && k > prev_k) {

				Expt[i].Pop[j].no_hosts += 1;
			}

			prev_i = i;
			prev_j = j;
			prev_k = k;
		}

		Data_File.close();

	/*allocate memory and sort data for each individual*/

		for (int i = 0; i < n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				Expt[i].Pop[j].Hosts = new host[Expt[i].Pop[j].no_hosts];
			}
		}


		Expt[0].Pop[0].Hosts[0].Data.no_measurement_times = 1;

		Data_File.open(filename);

		while (!Data_File.eof()) {

			getline(Data_File, experiment, '\t');
			getline(Data_File, group, '\t');
			getline(Data_File, individual, '\t');
			getline(Data_File, strain, '\t');
			getline(Data_File, refresh, '\t');
			getline(Data_File, dose, '\t');
			getline(Data_File, inoculum, '\t');
			getline(Data_File, time_measurement, '\t');
			getline(Data_File, time_infection, '\t');
			getline(Data_File, v, '\t');
			getline(Data_File, x, '\t');
			getline(Data_File, d, '\t');
			getline(Data_File, lod_v, '\t');
			getline(Data_File, lod_x, '\n');

			if (Data_File.eof()) break;

			i = stoi(experiment);
			j = stoi(group);
			k = stoi(individual);

			if (i != prev_i || j != prev_j || k != prev_k) {

				t = 0;

				Expt[i].Pop[j].Hosts[k].Data.no_measurement_times = 1;
			}

			if (t > 0 && i == prev_i && j == prev_j && k == prev_k) {

				Expt[i].Pop[j].Hosts[k].Data.no_measurement_times += 1;
			}


			Expt[i].Pop[j].Hosts[k].strain = strain;
			Expt[i].Pop[j].Hosts[k].refresh = refresh;
			Expt[i].Pop[j].Hosts[k].dose_size = stod(dose);
			Expt[i].Pop[j].Hosts[k].time_infection = stoi(time_infection);
			Expt[i].Pop[j].Hosts[k].viral_inoculum = stod(inoculum);
			Expt[i].Pop[j].Hosts[k].Data.lod_v = pow(10, stod(lod_v));/*limit of detection (extracellular virus)*/
			Expt[i].Pop[j].Hosts[k].Data.lod_x = pow(10, stod(lod_x));/*limit of detection (intracellular virus)*/

			Expt[i].Pop[j].Hosts[k].Data.Time[t] = stoi(time_measurement);

			if (v.compare("NA") == 0) {Expt[i].Pop[j].Hosts[k].Data.V[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.V[t] = stod(v);}

			if (x.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.X[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.X[t] = stod(x); }

			if (d.compare("NA") == 0) { Expt[i].Pop[j].Hosts[k].Data.D[t] = na; }
			else { Expt[i].Pop[j].Hosts[k].Data.D[t] = stod(d); }

			prev_i = i;
			prev_j = j;
			prev_k = k;

			t = t + 1;
		}


	Data_File.close();
}

void Read_Param_Values(string filename, model_type& model, fit_params& Fit, expt*Expt) {

	for (int i = 0; i < model.n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				ifstream Param_File;

				Param_File.open(filename);

				string param_name, sample, param_value_string;

				while (!Param_File.eof())
				{
					getline(Param_File, param_name, '\t');
					getline(Param_File, param_value_string, '\n');

						if (param_name == "n_expt")					model.n_expt = stoi(param_value_string);
	
						if (param_name == "beta")					Expt[i].Pop[j].Hosts[k].Params.beta = stod(param_value_string);
						if (param_name == "kappa_i")				Expt[i].Pop[j].Hosts[k].Params.kappa_i = stod(param_value_string);
						if (param_name == "kappa_t")				Expt[i].Pop[j].Hosts[k].Params.kappa_t = stod(param_value_string);
						if (param_name == "kappa_v")				Expt[i].Pop[j].Hosts[k].Params.kappa_v = stod(param_value_string);
						if (param_name == "omega")					Expt[i].Pop[j].Hosts[k].Params.omega = stod(param_value_string);
						if (param_name == "epsilon")				Expt[i].Pop[j].Hosts[k].Params.epsilon= stod(param_value_string);	
						if (param_name == "tau")					Expt[i].Pop[j].Hosts[k].Params.tau = stod(param_value_string);
						if (param_name == "s_n")					Expt[i].Pop[j].Hosts[k].Params.s_n = stod(param_value_string);
						if (param_name == "prop_v")					Expt[i].Pop[j].Hosts[k].Params.prop_v = stod(param_value_string);

						if (param_name == "in_vitro_model")			model.in_vitro_model = stoi(param_value_string);
						if (param_name == "drug_moa")				model.drug_moa = stoi(param_value_string);
						if (param_name == "hill_coeff")				Expt[i].Pop[j].Hosts[k].PK.hill_coeff = stod(param_value_string);
						if (param_name == "ec50")					Expt[i].Pop[j].Hosts[k].PK.ec50 = stod(param_value_string);
						if (param_name == "emax")					Expt[i].Pop[j].Hosts[k].PK.emax = stod(param_value_string);
								
						if (param_name == "no_fit")					Fit.no_fit = stoi(param_value_string);
						if (param_name == "no_particles")			Fit.no_particles = stoi(param_value_string);
						if (param_name == "iterations")				Fit.iterations = stod(param_value_string);
						if (param_name == "burn_in")				Fit.burn_in = stod(param_value_string);
						if (param_name == "thin")					Fit.thin = stoi(param_value_string);
						if (param_name == "scaling")				Fit.scaling = stod(param_value_string);
						if (param_name == "acpt_rate")				Fit.req_acpt_rate = stod(param_value_string);
						if (param_name == "no_post_samples")		Fit.no_post_samples = stoi(param_value_string);
						if (param_name == "block_update")			Fit.block_update = stoi(param_value_string);
						if (param_name == "cov_start")				Fit.cov_start = stoi(param_value_string);
						if (param_name == "cov_stop")				Fit.cov_stop = stoi(param_value_string);
						if (param_name == "cov_update")				Fit.cov_update = stoi(param_value_string);

						if (param_name == "error_v")				Expt[i].Pop[j].Hosts[k].Data.error_v = stod(param_value_string);
						if (param_name == "error_x")				Expt[i].Pop[j].Hosts[k].Data.error_x = stod(param_value_string);
						
				}
				
				Param_File.close();

			}
		}
	}
}

void Read_Param_Fit(string filename, int n_expt, fit_params& Fit, vector<string>& Param_Fit, vector<string>& Prior_Dist) {

	Fit.logscale = new int[Fit.no_fit];/*0-fitting parameter on linear scale, 1-fitting parameter on log10 scale*/
	Fit.level = new int[Fit.no_fit]; /*0-fitting at global level, 1- fitting at experiment level, 2-fitting at group level*/
	Fit.block = new int[Fit.no_fit];
	Fit.init_sd = new double[Fit.no_fit];/*initial jump size for proposing values*/
	Fit.max_sd = new double[Fit.no_fit];/*max jump size for proposing values*/
	Fit.min_sd = new double[Fit.no_fit];/*min jump size for proposing values*/
	Fit.start_values = new double[Fit.no_fit];/*starting value of each param we want to fit*/
	Fit.sample = new int[Fit.no_fit];

	ifstream Par_Fit;

	Par_Fit.open(filename);

	string hyperparameter,logscale, init_sd, sd_max, sd_min, level, block, param_name, param_value_string, distribution, lower, upper, sample;

	int i =0;
	int n = 0;

	vector<string> temp3;
	vector<string> temp4;

	while (!Par_Fit.eof())
	{
		
		if (n == Fit.no_fit) break;

		getline(Par_Fit, param_name, '\t');
		getline(Par_Fit, param_value_string, '\t');
		getline(Par_Fit, logscale, '\t');
		getline(Par_Fit, init_sd, '\t');
		getline(Par_Fit, sd_min, '\t');
		getline(Par_Fit, sd_max, '\t');
		getline(Par_Fit, level, '\t');
		getline(Par_Fit, block, '\t');
		getline(Par_Fit, distribution, '\t');
		getline(Par_Fit, lower, '\t');
		getline(Par_Fit, upper, '\n');

		//cout << param_name <<','<<param_value_string<< endl;

		if (i == 0) {

			Param_Fit.push_back(param_name);
			Prior_Dist.push_back(distribution);
		}

		Fit.start_values[n] = stod(param_value_string);
		Fit.logscale[n] = stoi(logscale);
		Fit.init_sd[n] = stod(init_sd);
		Fit.max_sd[n] = stod(sd_max);
		Fit.min_sd[n] = stod(sd_min);
		Fit.level[n] = stoi(level);
		Fit.block[n] = stoi(block);
			
		vector<string> temp1;
			
		vector<string> temp2;

			for (int j = 0; j < n_expt;j++) {

				temp1.push_back(lower);
				temp2.push_back(upper);
				
			}

			Fit.lower.push_back(temp1);
			Fit.upper.push_back(temp2);

		n = n + 1;
	}

	Par_Fit.close();
}

void Read_Dose_Schedule(string filename, expt& Expt) {

	ifstream Dose_Scheulde;
	Dose_Scheulde.open(filename);

	string hour, number;

	int no_doses;

	no_doses = -1;

	while (!Dose_Scheulde.eof()) {

		getline(Dose_Scheulde, number, '\t');
		getline(Dose_Scheulde, hour, '\n');

		no_doses += 1;

		if (Dose_Scheulde.eof()) break;
	}

	Dose_Scheulde.close();

	for (int j = 0; j < Expt.no_groups;j++) {

		for (int k = 0; k < Expt.Pop[j].no_hosts; k++) {

			Expt.Pop[j].Hosts[k].no_doses = no_doses;

			Expt.Pop[j].Hosts[k].Doses = new dose[no_doses];

		}
	}

	int t = 0;

	Dose_Scheulde.open(filename);

	while (!Dose_Scheulde.eof()) {

		getline(Dose_Scheulde, number, '\t');
		getline(Dose_Scheulde, hour, '\n');

		if (Dose_Scheulde.eof()) break;

		for (int j = 0; j < Expt.no_groups;j++) {

			for (int k = 0; k < Expt.Pop[j].no_hosts; k++) {

				Expt.Pop[j].Hosts[k].Doses[t].dose_no = stoi(number);
				Expt.Pop[j].Hosts[k].Doses[t].dose_time = stoi(hour);
				Expt.Pop[j].Hosts[k].Doses[t].dose_size = Expt.Pop[j].Hosts[k].dose_size;

				if (Expt.Pop[j].Hosts[k].Doses[t].dose_no == 0) { Expt.Pop[j].Hosts[k].first_dose = stoi(hour); }
				if (Expt.Pop[j].Hosts[k].Doses[t].dose_no == (Expt.Pop[j].Hosts[k].no_doses-1)) { Expt.Pop[j].Hosts[k].last_dose = stoi(hour); }

			}

		}

		t = t + 1;
	}

}

void Allocate_Memory(fit_params Fit, int n_expt, expt*Expt) {

	long length = long((Fit.iterations) / Fit.thin);

	long length_sub = long(((Fit.iterations - Fit.burn_in) / Fit.thin));

	for (int i = 0; i < n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int n = 0; n < Expt[i].Pop[j].no_hosts;n++) {

				Expt[i].Pop[j].Hosts[n].T = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].T[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) {
						Expt[i].Pop[j].Hosts[n].T[m][t] = new long double[RUNS];
					}
				}

				Expt[i].Pop[j].Hosts[n].E = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles; m++) {
					Expt[i].Pop[j].Hosts[n].E[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS; t++) { Expt[i].Pop[j].Hosts[n].E[m][t] = new long double[RUNS]; }
				}
				
				Expt[i].Pop[j].Hosts[n].I = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles; m++) {
					Expt[i].Pop[j].Hosts[n].I[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS; t++) { Expt[i].Pop[j].Hosts[n].I[m][t] = new long double[RUNS]; }
				}
				
				Expt[i].Pop[j].Hosts[n].X = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].X[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].X[m][t] = new long double[RUNS]; }
				}

				Expt[i].Pop[j].Hosts[n].V = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].V[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].V[m][t] = new long double[RUNS]; }
				}


				Expt[i].Pop[j].Hosts[n].D = new long double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].D[m] = new long double* [STEPS];
					for (int t = 0; t < STEPS;t++) { Expt[i].Pop[j].Hosts[n].D[m][t] = new long double[RUNS]; }
				}
	

				Expt[i].Pop[j].Hosts[n].LL = new double* [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) { Expt[i].Pop[j].Hosts[n].LL[m] = new double[DAYS]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.C = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.C[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_C = new double* [length_sub];
				for (int k = 0; k < length_sub;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_C[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.LC = new double[length];
				Expt[i].Pop[j].Hosts[n].Par_Chain.R0 = new double[length];

				Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_LC = new double[length_sub];
				Expt[i].Pop[j].Hosts[n].Par_Chain.Sub_R0 = new double[length_sub];

				Expt[i].Pop[j].Hosts[n].Par_Chain.SD = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.SD[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Posterior = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Posterior[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Prop_Acpt = new double* [length];
				for (int k = 0; k < length;k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Prop_Acpt[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Param_Sample = new double* [Fit.no_post_samples];
				for (int k = 0; k < Fit.no_post_samples; k++) { Expt[i].Pop[j].Hosts[n].Par_Chain.Param_Sample[k] = new double[Fit.no_fit]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Mean_Param_Sample = new double[Fit.no_fit];

				Expt[i].Pop[j].Hosts[n].Par_Chain.W = new double* [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles; m++) { Expt[i].Pop[j].Hosts[n].Par_Chain.W[m] = new double[DAYS]; }

				Expt[i].Pop[j].Hosts[n].Par_Chain.Lower = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.Median = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.Upper = new double[Fit.no_fit];

				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_val = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_val = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_sd = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_sd = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_posterior = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_posterior = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_ll = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_ll = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.curr_pr = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.new_pr = new double[Fit.no_fit];
				Expt[i].Pop[j].Hosts[n].Par_Chain.no_acpt = new double[Fit.no_fit];

				Expt[i].Pop[j].Hosts[n].Simulation = new host_predict[Fit.no_post_samples];

				for (int s = 0; s < Fit.no_post_samples;s++) {
					
					Expt[i].Pop[j].Hosts[n].Simulation[s].T = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].E = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].I = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].X = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].V = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].D = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].R = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Simulation[s].Eff = new double[STEPS];
				
				}

				Expt[i].Pop[j].Hosts[n].Bounds = new host_predict[3];

				for (int s = 0; s < 3;s++) {

					Expt[i].Pop[j].Hosts[n].Bounds[s].T = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].E = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].I = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].X = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].V = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].D = new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].R= new long double[STEPS];
					Expt[i].Pop[j].Hosts[n].Bounds[s].Eff = new double[STEPS];
					
				}

				
				Expt[i].Pop[j].Hosts[n].PK.A1 = new double** [Fit.no_particles];
				for (int m = 0; m < Fit.no_particles;m++) {
					Expt[i].Pop[j].Hosts[n].PK.A1[m] = new double* [STEPS];
					for (int t = 0; t < STEPS;t++) {
						Expt[i].Pop[j].Hosts[n].PK.A1[m][t] = new double[RUNS];
					}
				}
	
			}

			

		}
	
	}


}

void Init_Values(host& Host, int M, model_type model) {

	for (int m = 0; m < M;m++) {

		for (int r = 0; r < RUNS;r++) {
	
			Host.T[m][0][0] = 5000;
			Host.E[m][0][0] = 0;
			Host.I[m][0][0] = 0;
			Host.X[m][0][0] = 0;
			Host.D[m][0][0] = 0;
			
			Host.PK.A1[m][0][0] = 0;

			if (Host.time_infection==0 && RUN_IN==0) {

				Host.V[m][0][0] = Host.viral_inoculum;

			}

			else {

				Host.V[m][0][0] = 0;
			}

		}

		/*get carrying capacity*/

		Host.Params.K = (Host.Params.kappa_t * 14000) / ((Host.Params.s_n - Host.Params.kappa_t)); 

	}

}

void Drug_Dosage(host &Host, model_type model) {

	int t1, D;

	for (int t = 0; t < STEPS;t++) { Host.dosage[t] = 0; }

	D = Host.no_doses;/*total number of doses*/

	for (int d = 0; d < D;d++) {

		t1 = Host.Doses[d].dose_time+24*RUN_IN;

		Host.dosage[t1] += Host.Doses[d].dose_size;
	}
}

/*calculate mean across runs for each host, for each timestep*/
void Mean_Host(host_avg*Hosts_Avg, host*Host, int n, int runs, int steps) {

	for (int t = 0; t < steps; t++) {

		for (int r = 0; r < runs; r++) {

			Hosts_Avg[n].T[t][0] += Host[n].T[0][t][r] / double(runs);
			Hosts_Avg[n].E[t][0] += Host[n].E[0][t][r] / double(runs);
			Hosts_Avg[n].I[t][0] += Host[n].I[0][t][r] / double(runs);
			Hosts_Avg[n].X[t][0] += Host[n].X[0][t][r] / double(runs);
			Hosts_Avg[n].V[t][0] += Host[n].V[0][t][r] / double(runs);
			Hosts_Avg[n].D[t][0] += Host[n].D[0][t][r] / double(runs);

		}
	}
}

/*calculate standard deviation between runs for each timestep*/
void Std_Host(host_avg*Hosts_Avg, host*Host, int n, int runs, int steps) {

	double	diff_x,diff_v, diff_e, diff_i,diff_t,diff_d,
			var_x, var_v,var_d, var_e, var_i, var_t,
			sum_x, sum_v, sum_d, sum_e, sum_i, sum_t;

	for (int t = 0; t < steps;t++) {

		sum_x = 0; sum_t = 0; sum_v = 0; sum_d = 0; sum_e = 0; sum_i = 0; 

		for (int r = 0; r < runs; r++) {

			diff_t = Host[n].T[0][t][r] - Hosts_Avg[n].T[t][0];
			sum_t += pow(diff_t, 2);
			
			diff_e = Host[n].E[0][t][r] - Hosts_Avg[n].E[t][0];
			sum_e += pow(diff_e, 2);

			diff_i = Host[n].I[0][t][r] - Hosts_Avg[n].I[t][0];
			sum_i += pow(diff_i, 2);

			diff_x = Host[n].X[0][t][r] - Hosts_Avg[n].X[t][0];
			sum_x += pow(diff_x, 2); 
			
			diff_v = Host[n].V[0][t][r] - Hosts_Avg[n].V[t][0];
			sum_v += pow(diff_v, 2);

			diff_d = Host[n].D[0][t][r] - Hosts_Avg[n].D[t][0];
			sum_d += pow(diff_d, 2);

		}

		var_t = sum_t / double(runs);
		Hosts_Avg[n].T[t][1] = sqrt(var_t);

		var_e = sum_e / double(runs);
		Hosts_Avg[n].E[t][1] = sqrt(var_e);

		var_i = sum_i / double(runs);
		Hosts_Avg[n].I[t][1] = sqrt(var_i); 
		
		var_x = sum_x / double(runs);
		Hosts_Avg[n].X[t][1] = sqrt(var_x);

		var_v = sum_v / double(runs);
		Hosts_Avg[n].V[t][1] = sqrt(var_v);

		var_d = sum_d / double(runs);
		Hosts_Avg[n].D[t][1] = sqrt(var_d);

	}

}

void SubSample(double**C, double**S, double*LC, double*S2, double*R0, double*S3, int no_fit, long length, double burn_in, int thin) {

	long k2;

	int k0 = int(burn_in / thin);

	k2 = 0;

	for (int k = k0; k < length; k++) {

		S2[k2] = LC[k];
		S3[k2] = R0[k];

		for (int n = 0; n < no_fit; n++) {

			S[k2][n] = C[k][n];
		}

		k2 += 1;
	}
}

double Get_Quantile_Chains(double**C, int iterations, int n, double p) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	double*X;
	X = new double[iterations];

	for (int i = 0; i < iterations; i++) { X[i] = C[i][n]; }

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, p);

	return result;
}

double Get_Median(double*X, int iterations) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, 0.5);

	return result;
}

double Get_Lower(double*X, int iterations) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, 0.025);

	return result;
}

double Get_Upper(double*X, int iterations) {

	double temp, result;
	gsl_rng_set(g2, long(time(NULL)));

	/*sort values in ascending order*/

	for (int i = iterations - 1; i > 0; --i) {

		for (int j = 0; j < i; ++j) {

			if (X[j] > X[j + 1]) {

				temp = X[j];
				X[j] = X[j + 1];
				X[j + 1] = temp;

			}
		}
	}

	result = gsl_stats_quantile_from_sorted_data(X, 1, iterations, 0.975);

	return result;
}

void Set_New_Parameter_Values(host &Host, const vector<string> &v, int no_fit, double*S) {

	for (int n = 0; n < no_fit; n++) {
	
		if (strcmp(v[n].c_str(), "beta") == 0) { Host.Params.beta = S[n]; }
		if (strcmp(v[n].c_str(), "hill_coeff") == 0) { Host.PK.hill_coeff = S[n]; }
		if (strcmp(v[n].c_str(), "ec50") == 0) { Host.PK.ec50 = S[n]; }
		if (strcmp(v[n].c_str(), "omega") == 0) { Host.Params.omega = S[n]; }
		if (strcmp(v[n].c_str(), "prop_v") == 0) { Host.Params.prop_v = S[n]; }
		if (strcmp(v[n].c_str(), "error_v") == 0) { Host.Data.error_v= S[n]; }
		if (strcmp(v[n].c_str(), "error_x") == 0) { Host.Data.error_x = S[n]; }
	}
}

void Get_Quantile_Simulations(host& Host,host_predict*Simulation,host_predict Bound, int no_simulations, int T, double p){

	double temp;
	double*S;
	S = new double[no_simulations];

	for (int t = 0; t < T; t++) {

		/*sort values in ascending order*/

		for (int i = 0; i < no_simulations - 1; i++) {

			for (int j = i + 1; j < no_simulations; j++) {

				if (Simulation[i].T[t] > Simulation[j].T[t]) {

					temp = Simulation[i].T[t];
					Simulation[i].T[t] = Simulation[j].T[t];
					Simulation[j].T[t] = temp;

				}
			}
		}

		for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].T[t]; }

		/*get quantile at time t*/

		Bound.T[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);
		
		/*sort values in ascending order*/

		for (int i = 0; i < no_simulations - 1; i++) {

			for (int j = i + 1; j < no_simulations; j++) {

				if (Simulation[i].E[t] > Simulation[j].E[t]) {

					temp = Simulation[i].E[t];
					Simulation[i].E[t] = Simulation[j].E[t];
					Simulation[j].E[t] = temp;

				}
			}
		}

		for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].E[t];}

		/*get quantile at time t*/

		Bound.E[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

		/*sort values in ascending order*/

		for (int i = 0; i < no_simulations - 1; i++) {

			for (int j = i + 1; j < no_simulations; j++) {

				if (Simulation[i].I[t] > Simulation[j].I[t]) {

					temp = Simulation[i].I[t];
					Simulation[i].I[t] = Simulation[j].I[t];
					Simulation[j].I[t] = temp;

				}
			}
		}

		for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].I[t]; }

		/*get quantile at time t*/

		Bound.I[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);/*sort values in ascending order*/

		/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].X[t] > Simulation[j].X[t]) {

						temp = Simulation[i].X[t];
						Simulation[i].X[t] = Simulation[j].X[t];
						Simulation[j].X[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].X[t]; }

		/*get quatile at time t*/
			
			Bound.X[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

		/*sort values in ascending order*/
			
			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].V[t] > Simulation[j].V[t]) {

						temp = Simulation[i].V[t];
						Simulation[i].V[t] = Simulation[j].V[t];
						Simulation[j].V[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].V[t];}
			
		/*get quantile at time t*/
			
			Bound.V[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

		/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].D[t] > Simulation[j].D[t]) {

						temp = Simulation[i].D[t];
						Simulation[i].D[t] = Simulation[j].D[t];
						Simulation[j].D[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].D[t]; }

		/*get quantile at time t*/
			
			Bound.D[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

			/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].R[t] > Simulation[j].R[t]) {

						temp = Simulation[i].R[t];
						Simulation[i].R[t] = Simulation[j].R[t];
						Simulation[j].R[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].R[t]; }

			/*get quantile at time t*/

			Bound.R[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

			/*sort values in ascending order*/

			for (int i = 0; i < no_simulations - 1; i++) {

				for (int j = i + 1; j < no_simulations; j++) {

					if (Simulation[i].Eff[t] > Simulation[j].Eff[t]) {

						temp = Simulation[i].Eff[t];
						Simulation[i].Eff[t] = Simulation[j].Eff[t];
						Simulation[j].Eff[t] = temp;

					}
				}
			}

			for (int i = 0; i < no_simulations; i++) { S[i] = Simulation[i].Eff[t]; }

			/*get quantile at time t*/

			Bound.Eff[t] = gsl_stats_quantile_from_sorted_data(S, 1, no_simulations, p);

	}
}

void Find_Mean_Value(double**X, double *Y, int n, int m) {

	for (int i = 0; i < n; i++) {

		Y[i] = 0;

		for (int j = 0; j < m; j++) {

			//cout << X[j][i] << endl;

			Y[i] += X[j][i]/m;
		}

		//cout << Y[i] << endl;
	}

	//system("pause");

}

double Calculate_Rt(host& Host, int t, model_type model) {

	double Rt, effect, conc;
	double delta_prime, total_pop;
	int t1;

	t1 = t;

	conc = Host.D[0][t][0];

	if (conc > 0) { effect = Host.PK.emax / (1 + pow((Host.PK.ec50 / conc), Host.PK.hill_coeff)); }

	else { effect = 0; }

	total_pop = Host.T[0][t1][0] + Host.E[0][t1][0] + Host.I[0][t1][0];

	delta_prime = Host.Params.kappa_t * (1 + (total_pop/Host.Params.K)) +Host.Params.kappa_i;
	
	if (model.in_vitro_model == 0) {

		if (model.drug_moa == 1) {

				Rt = (Host.Params.beta * Host.T[0][t1][0] * Host.Params.omega * (1 - effect) * Host.Params.prop_v * Host.Params.tau) / (Host.Params.kappa_v * delta_prime * (Host.Params.tau * (1 - effect) + delta_prime));
			
			}

		else {
			
				Rt = (Host.Params.beta * Host.T[0][t1][0] * Host.Params.omega * (1 - effect) * Host.Params.prop_v * Host.Params.tau) / (Host.Params.kappa_v * delta_prime * (Host.Params.tau + delta_prime));
			
		}
	}

	else { Rt = (Host.Params.beta * Host.T[0][t1][0] * Host.Params.omega * (1 - effect) * Host.Params.prop_v) / (Host.Params.kappa_v * delta_prime); }
	
	//cout << (Host.Params.kappa_v * delta_prime * (Host.Params.tau * (1 - effect) + delta_prime)) << endl;

	//cout << Rt << endl;

	//system("pause");

	return Rt;

}

void Transform_Values(expt* Expt, fit_params Fit, int length, model_type model) {

	for (int i = 0; i < model.n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				for (int l = 0; l < length; l++) {

					for (int n = 0; n < Fit.no_fit; n++) {

						if (Fit.logscale[n] == 1) { Expt[i].Pop[j].Hosts[k].Par_Chain.C[l][n] = pow(10, Expt[i].Pop[j].Hosts[k].Par_Chain.C[l][n]); }

					}
				}
			}
		}
	}
}

void Calculate_Quantiles(expt* Expt, fit_params Fit, model_type model, int no_samples) {

	for (int i = 0; i < model.n_expt;i++) {

		for (int j = 0; j < Expt[i].no_groups;j++) {

			for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

				Expt[i].Pop[j].Hosts[k].Par_Chain.LL_Median = Get_Median(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, no_samples);
				Expt[i].Pop[j].Hosts[k].Par_Chain.R0_Median = Get_Median(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, no_samples);

				Expt[i].Pop[j].Hosts[k].Par_Chain.LL_Lower = Get_Lower(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, no_samples);
				Expt[i].Pop[j].Hosts[k].Par_Chain.R0_Lower = Get_Lower(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, no_samples);

				Expt[i].Pop[j].Hosts[k].Par_Chain.LL_Upper = Get_Upper(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, no_samples);
				Expt[i].Pop[j].Hosts[k].Par_Chain.R0_Upper = Get_Upper(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, no_samples);

				for (int n = 0; n < Fit.no_fit; n++) {

					Expt[i].Pop[j].Hosts[k].Par_Chain.Median[n] = Get_Quantile_Chains(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, no_samples, n, 0.5);
					Expt[i].Pop[j].Hosts[k].Par_Chain.Lower[n] = Get_Quantile_Chains(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, no_samples, n, 0.025);
					Expt[i].Pop[j].Hosts[k].Par_Chain.Upper[n] = Get_Quantile_Chains(Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, no_samples, n, 0.975);

				}
			}
		}
	}
}

void Posterior_Sample(expt*Expt, fit_params Fit, int iterations, model_type model) {

	int x;
	double* W;
	W = new double[iterations];

	for (int k = 0; k < iterations; k++) { W[k] = 1 / iterations; }
	
	/*all values equally likely to be sampled*/

	F1 = gsl_ran_discrete_preproc(iterations, W);

	for (int s = 0; s < Fit.no_post_samples; s++) {

		x = int(gsl_ran_discrete(g2, F1));

		for (int n = 0; n <Fit.no_fit; n++) {

			switch(Fit.level[n]){

			case 0:

				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[0].Pop[0].Hosts[0].Par_Chain.Sub_C[x][n];
							

						}
					}
				}

				
				break;

			case 1:

				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[i].Pop[0].Hosts[0].Par_Chain.Sub_C[x][n];
						}
					}
				}

				break;

			case 2:
		
				for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {


						for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

							Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s][n] = Expt[i].Pop[j].Hosts[0].Par_Chain.Sub_C[x][n];
						}
					}
				}

				break;

			}
		}
	}


	gsl_ran_discrete_free(F1);

	delete[] W;
}

double Calculate_DIC(expt*Expt,fit_params Fit, model_type model, const vector<string>&Param_Fit, param_vary& Param_Vary) {

	double pd, dic;

	pd = 0; dic = 0;

	/*For each estimated parameter, find mean value of posterior samples*/

		for (int i = 0; i < model.n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

					Find_Mean_Value(Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample, Expt[i].Pop[j].Hosts[k].Par_Chain.Mean_Param_Sample, Fit.no_fit, Fit.no_post_samples);
				}
			}
		}


	/*Calculate deviance at mean posterior values*/

		double dev_post_mean, dev_post_sample, mean_dev_post_sample, ll;

		dev_post_mean = 0;
		ll = 0;

		for (int i = 0; i < model.n_expt;i++) {

			for (int j = 0; j < Expt[i].no_groups;j++) {

				for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

					Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Mean_Param_Sample);
	
					ll=Exact_Likelihood(Fit, 0, 0, Expt[i].Pop[j].Hosts[k], Param_Fit, Param_Vary, model);
			
					dev_post_mean += -2 * ll;

				}
			}
		}


	/*Calculate mean of deviance evaluted at each individual posterior sample*/

		mean_dev_post_sample = 0;

		for (int s = 0; s < Fit.no_post_samples; s++) {

			dev_post_sample = 0;

			ll = 0;

			for (int i = 0; i < model.n_expt;i++) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

						Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s]);

						/*deterministic model - calculate loglikelihood exactly*/

						ll=Exact_Likelihood(Fit, 0, 0, Expt[i].Pop[j].Hosts[k], Param_Fit, Param_Vary, model);

						dev_post_sample += -2 * ll;
					}
				}
			}

			mean_dev_post_sample += dev_post_sample / Fit.no_post_samples;
		}

	/*Calculate DIC*/

		pd += mean_dev_post_sample - dev_post_mean;

		dic += mean_dev_post_sample + pd;

	return dic;
}

double Block_Size(expt* Expt, fit_params Fit, model_type model) {

	double size=0;

	for (int n = 0; n < Fit.no_fit; n++) {

		if (Fit.block[n] == 1) {

			if (Fit.level[n] == 0) { size += 1; }

			if (Fit.level[n] == 1) { size += 1*model.n_expt; }

			if (Fit.level[n] == 2) {

				for (int i = 0; i < model.n_expt;i++) {

					size += Expt[i].no_groups;
				}

			}
		}
	}

	return size;
}

void Covariance_Matrix(int size, int k1, int k2, expt* Expt, fit_params Fit, model_type model, gsl_matrix* Cov, gsl_matrix* Corr) {

	double corr, cov;

	double** temp; /*stores parameter values between iterations k1 and k2*/
	temp = new double* [size];

	int m = k2 - k1;

	for (int s = 0; s < size; s++) { temp[s] = new double[m]; }

	int s = 0;

	for (int n = 0; n < Fit.no_fit; n++) {

		switch (Fit.level[n]) {

		case 0:

			for (int k = k1; k < k2; k++) {
	
				temp[s][k - k1] = Expt[0].Pop[0].Hosts[0].Par_Chain.C[k][n];
	
				
			}

			s = s + 1;

			break;

		case 1:

			for (int i = 0; i < model.n_expt;i++) {

				for (int k = k1; k < k2; k++) {

					temp[s][k - k1] = Expt[i].Pop[0].Hosts[0].Par_Chain.C[k][n];
				}

				s = s + 1;
			}

			break;

		case 2:

			for (int i = 0; i < model.n_expt;i++) {

				for (int j = 0; j < Expt[i].no_groups;j++) {

					for (int k = k1; k < k2; k++) {

						temp[s][k - k1] = Expt[i].Pop[j].Hosts[0].Par_Chain.C[k][n];
					}

					s = s + 1;
				}
			}

			break;

		}

	}

	/*compute covariance/correlation matrix*/

	for (int s1 = 0; s1 < size; s1++) {

		for (int s2 = 0; s2 < size; s2++) {

			gsl_matrix_set(Corr, s1, s2, 0); /*reset correlation matrix to zero*/
			corr = gsl_stats_correlation(temp[s1], 1, temp[s2], 1, m); /*computes correlation efficient between temp[n] and temp[m]*/
			gsl_matrix_set(Corr, s1, s2, corr); /*populate correlation matrix*/

			gsl_matrix_set(Cov, s1, s2, 0);/*reset covariance matrix to zero*/
			cov = gsl_stats_covariance(temp[s1], 1, temp[s2], 1, m);/*computes covariance between temp[n] and temp[m]*/
			gsl_matrix_set(Cov, s1, s2, cov);/*populate covariance matrix*/

		}
	}

	for (int s = 0; s < size; s++) { delete temp[s]; }

	delete[] temp;
}

