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
#include<limits.h>
#include <omp.h>
#include"randlib_par.h"
#include"Structures.h"
#include"Functions.h"
#include"In_Vitro_Model.h"
#include"MCMC.h"

using namespace std;

void In_Vitro_Model(host& Host, int thread, int start, int steps, int runs, model_type model) {

	double
		h_t, h_e, h_i, h_v, h_x,
		p_t, p_e, p_i, p_v, p_x,
		t1;

	long double
		current_t, current_e, current_i, current_x, current_v, current_d, total_cells,
		new_e, new_x, new_i, new_v, innoc_v,
		rep_t,
		out_t, out_e, out_x, out_i,
		deaths_t, deaths_e, deaths_y, extra_deaths_y, deaths_i, deaths_v;

	int m = thread;

	model_params params = Host.Params;

	for (int r = 0; r < runs; r++) {

		for (int t = start; t < steps + 1; t++)

		{

			t1 = int(t / (1 / DT)); /*day number*/

			/*current population in each compartment*/

			current_t = Host.T[m][t][r]; /*target cells*/
			current_e = Host.E[m][t][r]; /*exposed cells*/
			current_i = Host.I[m][t][r]; /*infectious cells*/
			current_x = Host.X[m][t][r]; /*intracellular virus*/
			current_v = Host.V[m][t][r]; /*extracellular virus*/
			current_d = Host.D[m][t][r]; /*drug*/

			total_cells = current_t + current_e + current_i;

			/*hazard of leaving each compartment*/

			h_t = params.kappa_t * (1 + (total_cells / params.K)) + Infection_Rate(Host, thread, t, r, model, params.beta, Host.PK.emax, Host.PK.ec50, Host.PK.hill_coeff);
			h_e = params.kappa_t * (1 + (total_cells / params.K)) + params.kappa_i + Progression_Rate(Host, thread, t, r, model, params.tau, Host.PK.emax, Host.PK.ec50, Host.PK.hill_coeff);
			h_i = params.kappa_t * (1 + (total_cells / params.K)) + params.kappa_i;
			h_x = params.epsilon;
			h_v = params.kappa_v;

			/*convert hazards to probabilities*/

			p_t = 1 - exp(-h_t * DT);
			p_e = 1 - exp(-h_e * DT);
			p_i = 1 - exp(-h_i * DT);
			p_x = 1 - exp(-h_x * DT);
			p_v = 1 - exp(-h_v * DT);

			/*determine movement between compartments*/

				rep_t = params.s_n * current_t * DT;
				out_t = current_t * p_t;
				if (h_t > 0) { deaths_t = out_t * ((params.kappa_t * (1 + (total_cells / params.K))) / (h_t)); }
				else { deaths_t = 0; }

				if (model.in_vitro_model == 0) {

					/*with latent period*/

					new_e = out_t - deaths_t;
					out_e = current_e * p_e;
					deaths_e = out_e * ((params.kappa_t * (1 + (total_cells / params.K)) + params.kappa_i) / h_e);

					new_i = out_e - deaths_e;
				}

				else {

					/*without latent period*/

					new_e = 0;
					out_e = 0;
					deaths_e = 0;

					new_i = out_t - deaths_t;

				}

				out_i = current_i * p_i;
				
				new_x = Virus_Proliferation(Host, thread, t, r, model, params.omega, Host.PK.emax, Host.PK.ec50, Host.PK.hill_coeff);
				out_x = current_x * p_x;

				new_v = params.prop_v * out_x;
				deaths_v = current_v * p_v;
	
				if (t == (Host.time_infection + 24 * RUN_IN - 1) && t > 0) { innoc_v = Host.viral_inoculum; }
				else { innoc_v = 0; }

			/*calculate population size at next timestep*/

				if (model.in_vitro_model == 0) {

					/*with latent period*/

					Host.T[m][t + 1][r] = current_t + rep_t - deaths_t - new_e;
					Host.E[m][t + 1][r] = current_e + new_e - deaths_e - new_i;
					Host.I[m][t + 1][r] = current_i + new_i - out_i;
					Host.X[m][t + 1][r] = current_x + new_x - out_x;
					Host.V[m][t + 1][r] = current_v + innoc_v + new_v - deaths_v;
					Host.D[m][t + 1][r] = PK_Model(Host, m, t, r, model);
				}

				else {

					/*without latent period*/

					Host.T[m][t + 1][r] = current_t + rep_t - deaths_t - new_i;
					Host.E[m][t + 1][r] = 0;
					Host.I[m][t + 1][r] = current_i + new_i - out_i;
					Host.X[m][t + 1][r] = current_x + new_x - out_x;
					Host.V[m][t + 1][r] = current_v + innoc_v + new_v - deaths_v;
					Host.D[m][t + 1][r] = PK_Model(Host, m, t, r, model);

				}
		}

	}

}

double Infection_Rate(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff) {

	double new_rate, conc, effect;
	int m = thread;

	model_params params = Host.Params;

	conc = Host.D[m][t][r];

	if (conc > 0) {

		if (model.drug_moa == 0) { effect = emax / (1 + pow((ec50 / conc), hill_coeff)); }

		else { effect = 0; }
	}

	else {

		effect = 0;
	}

	new_rate = Host.V[m][t][r] * rate * (1 - effect);

	return new_rate;
}

double Progression_Rate(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff) {

	double new_rate, conc, effect;
	int m = thread;

	model_params params = Host.Params;

	conc = Host.D[m][t][r];

	if (conc > 0) {

		if (model.drug_moa == 1) { effect = emax / (1 + pow((ec50 / conc), hill_coeff)); }

		else { effect = 0; }
	}

	else {

		effect = 0;
	}

	new_rate = rate * (1 - effect);

	return new_rate;
}

long double Virus_Proliferation(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff) {

	long double new_v, conc, effect;
	int m = thread;
	int t1;

	model_params params = Host.Params;

	t1 = t;

	conc = Host.D[m][t][r];

	if (conc > 0) {

		if (model.drug_moa == 2) { effect = emax / (1 + pow((ec50 / conc), hill_coeff)); }

		else { effect = 0; }
	}

	else {

		effect = 0;
	}

	new_v = rate * (1 - effect) * Host.I[m][t1][r] * DT;
		
	return new_v;
}

double PK_Model(host &Host, int thread, int t, int r, model_type model) {

	int m = thread;

	Host.PK.A1[m][t + 1][r] = Host.PK.A1[m][t][r] + Host.dosage[t];

	return Host.PK.A1[m][t + 1][r];

}
