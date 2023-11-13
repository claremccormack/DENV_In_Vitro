//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																																									//
// In Vitro Model (Modelling the impact of JNJ-1802, a first-in-class dengue inhibitor blocking the NS3-NS4B interaction, on in-vitro DENV-2 dynamics)				//
//																																									//
//																																									//
// The model describes the dynamics of infection in an individual well. 																							//
// Each individual well is part of a larger group of wells receiving a particular concentration of JNJ-1802. 														//							
// The model is fitted to data generated from multiple experiments, and thus the model is structured as follows:													//
//																																									//
// Experiment -> Treatment Group -> Individual Well																													//
//																																									//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
#include<omp.h>
#include"randlib_par.h"
#include"Structures.h"
#include"In_Vitro_Model.h"
#include"Functions.h"
#include"MCMC.h"
#include"gsl/gsl_randist.h"
#include"gsl/gsl_rng.h"
#include"gsl/gsl_matrix.h"
#include"gsl/gsl_linalg.h"
#include"gsl/gsl_cdf.h"
#include"gsl/gsl_math.h"

using namespace std;

const double na = sqrt(-1);

int main(int argc, char*argv[]) {

	//////////////////
	//*Model Setup*///
	//////////////////

		fit_params Fit;
		param_vary Param_Vary;
		model_type model;

			/*read in data from file*/

				expt* Expt;
				
				Read_Data(argv[1], Expt, model.n_expt);

			/*read in fixed parameter values from file*/

				Read_Param_Values(argv[2], model, Fit, Expt);

				Allocate_Memory(Fit, model.n_expt, Expt);

			/*read in list of parameter values we want to fit*/
		
				vector<string>Param_Fit;/*stores names of parameters we want to estimate*/
				vector<string>Prior_Dist; /*stores prior distribution type of parameters we want to estimate*/
		
				Read_Param_Fit(argv[3], model.n_expt, Fit, Param_Fit, Prior_Dist);

			/*read in dosing scheule for each experiment*/

				int n = 4;

				for (int i = 0; i < model.n_expt;i++) {

					Read_Dose_Schedule(argv[n], Expt[i]);

					n = n + 1;
				}


			/*read in initial covariance matrix*/

				gsl_matrix* Init_Cov;

				double size = Block_Size(Expt, Fit, model);

				Init_Cov = gsl_matrix_alloc(size,size);
			
				gsl_matrix_set_identity(Init_Cov);

		/*set seeds for random number generation*/

			initSeeds(unsigned(time(NULL)), unsigned(time(NULL)));

			srand(unsigned(time(NULL)));

			random_device rd;
			mt19937 gen(rd());
	
		/*output data to file*/

					/*create directory for output files*/

					stringstream ss_m; /*extension for filenames*/

					ss_m << model.n_expt << ',' << Fit.no_fit << ',' << Param_Fit[0] << ',' << Fit.level[0] << ',' << Fit.start_values[0] << ',' << Fit.init_sd[0];

					_mkdir("Results"); /*make results directory*/

					ofstream results_data;
					results_data.open("Results/Data " + ss_m.str() + ".csv");

					results_data << "Experiment"<<','<<"Group" << ',' << "Well" << ',' << "Strain" << ',' << "Refresh"<<','<<"Concentration" << ',' << "Inoculum" << ',' << "Time" << ',' << "Virus_Extra" << ',' << "Virus_Intra"<<endl;

					for (int i = 0; i < model.n_expt;i++) {

						for (int j = 0; j < Expt[i].no_groups;j++) {

							for (int n = 0; n < Expt[i].Pop[j].no_hosts; n++) {

								for (int t = 0; t < Expt[i].Pop[j].Hosts[n].Data.no_measurement_times;t++) {

									results_data << i << ',' <<j<<','<< n << ',' << Expt[i].Pop[j].Hosts[n].strain << ',' << Expt[i].Pop[j].Hosts[n].refresh <<','<<Expt[i].Pop[j].Hosts[n].dose_size << ',' << Expt[i].Pop[j].Hosts[n].viral_inoculum << ',' << Expt[i].Pop[j].Hosts[n].Data.Time[t] << ',';

									if (isnan(Expt[i].Pop[j].Hosts[n].Data.V[t]) == 1) {

										results_data << "NA" << ',';

									}

									else {

										results_data << Expt[i].Pop[j].Hosts[n].Data.V[t] << ',';
									}

									if (isnan(Expt[i].Pop[j].Hosts[n].Data.X[t]) == 1) {

										results_data << "NA" << endl;

									}

									else {

										results_data << Expt[i].Pop[j].Hosts[n].Data.X[t] << endl;
									}

								}

								results_data << endl;
							}

							results_data << endl;
						}
				}


	//////////////////////////
	///*MODEL FITTING*///////
	/////////////////////////

				/*populate arrays with initial values*/

					for (int i = 0; i < model.n_expt;i++) {

						for (int j = 0; j < Expt[i].no_groups;j++) {

							for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {
							
								for (int n = 0; n < Fit.no_fit; n++) {

									Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] = Fit.start_values[n];							
									Expt[i].Pop[j].Hosts[k].Par_Chain.curr_val[n] = Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n];
	
									Expt[i].Pop[j].Hosts[k].Par_Chain.SD[0][n] = Fit.init_sd[n];
									Expt[i].Pop[j].Hosts[k].Par_Chain.curr_sd[n] = Fit.init_sd[n];

									if (Fit.logscale[n] == 1) {

										Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n] = log10(Expt[i].Pop[j].Hosts[k].Par_Chain.C[0][n]);
										Expt[i].Pop[j].Hosts[k].Par_Chain.curr_val[n] = log10(Expt[i].Pop[j].Hosts[k].Par_Chain.curr_val[n]);


									}
				
								}
							}
						}
					}

				/*run MCMC algorithm*/

					int length = long((Fit.iterations) / Fit.thin); /*Thin MCMC chains by factor Fit.thin*/

					int no_samples = int((Fit.iterations - Fit.burn_in) / Fit.thin);

					Run_MCMC(Fit, model.n_expt, Expt, Param_Fit, Prior_Dist, Param_Vary, model, Init_Cov);

					/*if working on log10 scale, transform values back to linear scale*/
						
					Transform_Values(Expt, Fit, length, model);
				
					for (int i = 0; i < model.n_expt; i++) {

							for (int j = 0; j < Expt[i].no_groups; j++) {

								for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

									SubSample(Expt[i].Pop[j].Hosts[k].Par_Chain.C, Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_C, Expt[i].Pop[j].Hosts[k].Par_Chain.LC, Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_LC, Expt[i].Pop[j].Hosts[k].Par_Chain.R0, Expt[i].Pop[j].Hosts[k].Par_Chain.Sub_R0, Fit.no_fit, length, Fit.burn_in, Fit.thin);

								}
							}
						}
			
				/*calculate quantiles for posterior estimates of parameter values*/

					Calculate_Quantiles(Expt, Fit, model, no_samples);

	///////////////////////////////////////////
	/////*POSTERIOR PREDICTIVE ANALYSIS*///////
	//////////////////////////////////////////

				/*Sample values from posterior distributions, accounting for level the parameter is fitted at*/

					Posterior_Sample(Expt, Fit, no_samples, model);

				/*For each posterior sample, run model once and store results*/

					for (int i = 0; i < model.n_expt;i++) {

						for (int j = 0; j < Expt[i].no_groups;j++) {

							for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

								for (int s = 0; s < Fit.no_post_samples; s++) {

									/*set parameter values*/

									Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s]);

									/*run model*/
								
										Drug_Dosage(Expt[i].Pop[j].Hosts[k], model);

										Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Param_Sample[s]);

										Init_Values(Expt[i].Pop[j].Hosts[k], Fit.no_particles, model);

										In_Vitro_Model(Expt[i].Pop[j].Hosts[k], 0, 0, STEPS - 2, 1, model);

										/*store results*/

										for (int t = 0; t < STEPS; t++) {

											Expt[i].Pop[j].Hosts[k].Simulation[s].T[t] = Expt[i].Pop[j].Hosts[k].T[0][t][0];
											Expt[i].Pop[j].Hosts[k].Simulation[s].E[t] = Expt[i].Pop[j].Hosts[k].E[0][t][0];
											Expt[i].Pop[j].Hosts[k].Simulation[s].I[t] = Expt[i].Pop[j].Hosts[k].I[0][t][0];
											Expt[i].Pop[j].Hosts[k].Simulation[s].X[t] = Expt[i].Pop[j].Hosts[k].X[0][t][0];
											Expt[i].Pop[j].Hosts[k].Simulation[s].V[t] = Expt[i].Pop[j].Hosts[k].V[0][t][0];
											Expt[i].Pop[j].Hosts[k].Simulation[s].D[t] = Expt[i].Pop[j].Hosts[k].D[0][t][0];
											Expt[i].Pop[j].Hosts[k].Simulation[s].R[t] = Calculate_Rt(Expt[i].Pop[j].Hosts[k], t, model);
											if (Expt[i].Pop[j].Hosts[k].Simulation[s].D[t] > 0) { Expt[i].Pop[j].Hosts[k].Simulation[s].Eff[t] = Expt[i].Pop[j].Hosts[k].PK.emax / (1 + pow((Expt[i].Pop[j].Hosts[k].PK.ec50 / Expt[i].Pop[j].Hosts[k].Simulation[s].D[t]), Expt[i].Pop[j].Hosts[k].PK.hill_coeff)); }
											else { Expt[i].Pop[j].Hosts[k].Simulation[s].Eff[t] = 0; }
										}
								}
								
								/*calculate quantiles across simulations*/

								Get_Quantile_Simulations(Expt[i].Pop[j].Hosts[k],Expt[i].Pop[j].Hosts[k].Simulation, Expt[i].Pop[j].Hosts[k].Bounds[0], Fit.no_post_samples, STEPS, .5);
								Get_Quantile_Simulations(Expt[i].Pop[j].Hosts[k],Expt[i].Pop[j].Hosts[k].Simulation, Expt[i].Pop[j].Hosts[k].Bounds[1], Fit.no_post_samples, STEPS, .025);
								Get_Quantile_Simulations(Expt[i].Pop[j].Hosts[k],Expt[i].Pop[j].Hosts[k].Simulation, Expt[i].Pop[j].Hosts[k].Bounds[2], Fit.no_post_samples, STEPS, .975);

						}
							
					}
				}
		

		////////////////////////////
		/////*Summary Statistics*///
		///////////////////////////

				/*run model with posterior median parameter values*/

					for (int i = 0; i < model.n_expt;i++) {

						for (int j = 0; j < Expt[i].no_groups;j++) {

							for (int k = 0; k < Expt[i].Pop[j].no_hosts; k++) {

								/*set parameter values*/

								Set_New_Parameter_Values(Expt[i].Pop[j].Hosts[k], Param_Fit, Fit.no_fit, Expt[i].Pop[j].Hosts[k].Par_Chain.Median);

								Expt[i].Pop[j].Hosts[k].Params.R0 = Calculate_Rt(Expt[i].Pop[j].Hosts[k],0,model);

								/*run model*/

								Init_Values(Expt[i].Pop[j].Hosts[k], Fit.no_particles, model);

								Drug_Dosage(Expt[i].Pop[j].Hosts[k],model);

								In_Vitro_Model(Expt[i].Pop[j].Hosts[k], 0, 0, STEPS - 2, RUNS, model);
							}
						}
					}

					ofstream results1;
					results1.open("Results/Hosts " + ss_m.str() + ".csv");

					for (int i = 0; i < model.n_expt; i++) {

						for (int j = 0; j < Expt[i].no_groups; j++) {

							for (int n = 0; n < Expt[i].Pop[j].no_hosts; n++) {

								for (int r = 0; r < 1; r++) {

									for (int t = 0; t < STEPS; t++) {

										results1 << Expt[i].Pop[j].Hosts[n].T[0][t][r] << ',' << Expt[i].Pop[j].Hosts[n].E[0][t][r] << ',' << Expt[i].Pop[j].Hosts[n].I[0][t][r] << ',' << Expt[i].Pop[j].Hosts[n].X[0][t][r] << ',' << Expt[i].Pop[j].Hosts[n].V[0][t][r] << ',' << Expt[i].Pop[j].Hosts[n].D[0][t][r] << endl; 			
									}

									results1 << endl;
								}

								results1 << endl;
							}

							results1 << endl;
						}

						results1 << endl;
					}

	//////////////////////////////////////////////////////////
	///*Calculate Deviance Information Criterion (DIC)*///////
	//////////////////////////////////////////////////////////

					double dic;

					dic = Calculate_DIC(Expt, Fit, model, Param_Fit, Param_Vary);

	//////////////////////
	///*Output Results*///
	//////////////////////
				
				ofstream results_mcmc;
				results_mcmc.open("Results/MCMC " + ss_m.str() + ".csv");

				for (int k = 0; k < length; k++) {

					for (int n = 0; n < Fit.no_fit; n++) {

						if (Fit.level[n] == 0) {

							if (Fit.logscale[n] == 1) {

								results_mcmc << log10(Expt[0].Pop[0].Hosts[0].Par_Chain.C[k][n]) << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.LC[k] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.R0[k] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.SD[k][n] << ',';

							}

							else {

								results_mcmc << Expt[0].Pop[0].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.LC[k] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.R0[k] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.SD[k][n] << ',';

							}

						}

						if (Fit.level[n] == 1) {

							for (int i = 0; i < model.n_expt;i++) {

									if (Fit.logscale[n] == 1) {

										results_mcmc << log10(Expt[i].Pop[0].Hosts[0].Par_Chain.C[k][n]) << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.LC[k] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.R0[k] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.SD[k][n] << ',';

									}

									else {

										results_mcmc << Expt[i].Pop[0].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.LC[k] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.R0[k] << ',' << Expt[i].Pop[0].Hosts[0].Par_Chain.SD[k][n] << ',';

									}
								}
						}

						if (Fit.level[n] == 2) {

							for (int i = 0; i < model.n_expt;i++) {

								for (int j = 0; j < Expt[i].no_groups;j++) {

										if (Fit.logscale[n] == 1) {

											results_mcmc << log10(Expt[i].Pop[j].Hosts[0].Par_Chain.C[k][n]) << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.LC[k] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.R0[k] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.SD[k][n] << ',';

										}

										else {

											results_mcmc << Expt[i].Pop[j].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.C[k][n] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.LC[k] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.R0[k] << ',' << Expt[i].Pop[j].Hosts[0].Par_Chain.SD[k][n] << ',';

										}
									}
							}
						}

						if (Fit.level[n] == 3) {

							for (int i = 0; i < model.n_expt;i++) {

								for (int j = 0; j < Expt[i].no_groups;j++) {

									for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

										if (Fit.logscale[n] == 1) {

											results_mcmc << log10(Expt[i].Pop[j].Hosts[l].Par_Chain.C[k][n]) << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.C[k][n] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.LC[k] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.R0[k] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.SD[k][n] << ',';

										}

										else {

											results_mcmc << Expt[i].Pop[j].Hosts[l].Par_Chain.C[k][n] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.C[k][n] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.LC[k] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.R0[k] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.SD[k][n] << ',';

										}
									}
								}
							}
						}
					}

					results_mcmc << endl;
				}

				ofstream results_mcmc_sub;
				results_mcmc_sub.open("Results/MCMC Sub " + ss_m.str() + ".csv");

					int k0 = int((Fit.iterations - Fit.burn_in) / Fit.thin);

					for (int k = 0; k < k0; k++) {

						for (int n = 0; n < Fit.no_fit; n++) {

							if (Fit.level[n] == 0) { results_mcmc_sub << Expt[0].Pop[0].Hosts[0].Par_Chain.Sub_C[k][n] << ','; }

							if (Fit.level[n] == 1) {

								for (int i = 0; i < model.n_expt;i++) {

										results_mcmc_sub << Expt[i].Pop[0].Hosts[0].Par_Chain.Sub_C[k][n] << ',';

								}
							}

							if (Fit.level[n] == 2) {

								for (int i = 0; i < model.n_expt; i++) {

									for (int j = 0; j < Expt[i].no_groups;j++) {

										results_mcmc_sub << Expt[i].Pop[j].Hosts[0].Par_Chain.Sub_C[k][n] << ',';
									}
								}
							}

							if (Fit.level[n] == 3) {

								for (int i = 0; i < model.n_expt;i++) {

									for (int j = 0; j < Expt[i].no_groups;j++) {

										for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

											results_mcmc_sub << Expt[i].Pop[j].Hosts[l].Par_Chain.Sub_C[k][n] << ',';
										}
									}
								}
							}
						}

						results_mcmc_sub << Expt[0].Pop[0].Hosts[0].Par_Chain.SD[k0-1][0]<<','<<Expt[0].Pop[0].Hosts[0].Par_Chain.Sub_LC[k] << ',' << Expt[0].Pop[0].Hosts[0].Par_Chain.Sub_R0[k] << endl;
					}

					results_mcmc_sub << endl;

				ofstream results_fit;
				results_fit.open("Results/MCMC Fit " + ss_m.str() + ".csv");

					k0 = int((Fit.iterations) / Fit.thin);

					for (int i = 0; i < model.n_expt;i++) {

						for (int j = 0; j < Expt[i].no_groups;j++) {

							for (int l = 0; l < Expt[i].Pop[j].no_hosts; l++) {

								for (int n = 0; n < Fit.no_fit; n++) {

									results_fit << Fit.burn_in << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.Prop_Acpt[k0 - 1][n] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.Median[n] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.Lower[n] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.Upper[n] << ',' << dic << endl;
								}

								results_fit << Fit.burn_in << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.Prop_Acpt[k0 - 1][0] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.LL_Median << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.LL_Lower << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.LL_Upper << ','<< dic << endl;
								results_fit << Fit.burn_in << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.Prop_Acpt[k0 - 1][0] << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.R0_Median << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.R0_Lower << ',' << Expt[i].Pop[j].Hosts[l].Par_Chain.R0_Upper << ','<< dic << endl;

								results_fit << endl;
							}

							results_fit << endl;
						}

						results_fit << endl;
					}

				ofstream results2;
				results2.open("Results/Bounds " + ss_m.str() + ".csv");

					for (int i = 0; i < model.n_expt;i++) {

					for (int j = 0; j < Expt[i].no_groups;j++) {

						for (int n = 0; n < Expt[i].Pop[j].no_hosts;n++) {

							for (int t = 0; t < STEPS; t++) {

								results2
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].T[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].T[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].T[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].E[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].E[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].E[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].I[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].I[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].I[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].X[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].X[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].X[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].V[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].V[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].V[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].D[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].D[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].D[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].R[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].R[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].R[t] << ','
									<< Expt[i].Pop[j].Hosts[n].Bounds[0].Eff[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[1].Eff[t] << ',' << Expt[i].Pop[j].Hosts[n].Bounds[2].Eff[t] << endl;

							}

							results2 << endl;
						}

						results2 << endl;
					}

					results2 << endl;
				}

}