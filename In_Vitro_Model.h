#ifndef IN_VITRO_MODEL_H
#define IN_VITRO_MODE_H

#include "Structures.h"

void In_Vitro_Model(host& Host, int thread, int start, int steps, int runs, model_type model);
double Infection_Rate(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff);
double Progression_Rate(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff);
long double Virus_Proliferation(host Host, int thread, int t, int r, model_type model, double rate, double emax, double ec50, double hill_coeff);
double PK_Model(host& Host, int thread, int t, int r, model_type model);

#endif
