#include <stdio.h>
#include <fstream>


#include "Strucrs.h"

int main(int argc, char *argv[]){

	unsigned int k;

	std::string buffer;
	std::ifstream IN_file;

	FILE *rhoF, *outF;

	Parameters param;

	param.IN_filename = argv[1];
	param.OUT_filename = argv[2];
	param.RHO_filename = argv[3];

	IN_file.open(param.IN_filename);

	IN_file >> param.N;
	std::getline(IN_file, buffer);

	IN_file >> param.dTau;
	std::getline(IN_file, buffer);

	IN_file >> param.end;
	std::getline(IN_file, buffer);

	IN_file >> param.n;
	std::getline(IN_file, buffer);

	IN_file >> param.omega;
	std::getline(IN_file, buffer);

	IN_file >> param.kappa;
	std::getline(IN_file, buffer);

	IN_file.close();

	param.omega *= pow(PI, 2) / 2;

	double *x = new double[param.N+1];
	double *psi_r = new double[param.N+1];
	double *psi_i = new double[param.N+1];
	double *H_r = new double[param.N+1];
	double *H_i = new double[param.N+1];

	double energia = 0.0, norma = 0.0, x_mean = 0.0, licznik=0.0;

	param.tau = 0.0;

	param.dx = 1.0 / param.N;

	for (k = 0; k < param.N+1; k++){

		x[k] = k*param.dx;

		psi_r[k] = sqrt(2.0)*sin(param.n*PI*x[k]);
		psi_i[k] = 0.0;

	}

	for (k = 1; k < param.N; k++) {

		H_r[k] = -(psi_r[k + 1] + psi_r[k - 1] - 2.0*psi_r[k]) / (2.0*pow(param.dx, 2)) + param.kappa*(x[k] - 0.5)*psi_r[k] * sin(param.omega*param.tau);
		H_i[k] = -(psi_i[k + 1] + psi_i[k - 1] - 2.0*psi_i[k]) / (2.0*pow(param.dx, 2)) + param.kappa*(x[k] - 0.5)*psi_i[k] * sin(param.omega*param.tau);
	}

	H_r[0] = 0.0;
	H_r[param.N] = 0.0;

	H_i[0] = 0.0;
	H_i[param.N] = 0.0;

	fopen_s(&rhoF, param.RHO_filename.c_str(), "wb");
	fopen_s(&outF, param.OUT_filename.c_str(), "wb");
	fprintf(outF, "t\tN\tx\tE\tsin(t)\n");

for (param.tau = 0.0; param.tau < param.end; param.tau += param.dTau) {

		norma = 0;
		energia = 0;
		x_mean = 0;
		if (licznik > param.end / 100) {
			printf("\r(%.2f/%.2f)", param.tau, param.end);

			for (k = 0; k < param.N + 1; k++) {
				fprintf(rhoF, "%e\t", pow(psi_r[k], 2) + pow(psi_i[k], 2));
				norma += pow(psi_r[k], 2) + pow(psi_i[k], 2);
				x_mean += x[k] * (pow(psi_r[k], 2) + pow(psi_i[k], 2));
				energia += H_r[k] * psi_r[k] + H_i[k] * psi_i[k];
			}

			fprintf(rhoF, "\n");
			fprintf(outF, "%e\t%e\t%e\t%e\t%e\n", param.tau, norma*param.dx, x_mean*param.dx, energia*param.dx, param.kappa*sin(param.omega*(param.tau)));

			licznik = 0.0;
		}

		 
		for (k = 0; k < param.N+1; k++)
			psi_r[k] = psi_r[k] + H_i[k] * param.dTau / 2;

		for (k = 1; k < param.N; k++)
			H_r[k] = -(psi_r[k + 1] + psi_r[k - 1] - 2.0*psi_r[k]) / (2.0*pow(param.dx, 2))
				+ param.kappa*(x[k] - 0.5)*psi_r[k] * sin(param.omega*(param.tau+param.dTau/2));

		for (k = 0; k < param.N+1; k++)
			psi_i[k] = psi_i[k] - H_r[k] * param.dTau;

		for (k = 1; k < param.N; k++)
			H_i[k] = -(psi_i[k + 1] + psi_i[k - 1] - 2.0*psi_i[k]) / (2.0*pow(param.dx, 2))
				+ param.kappa*(x[k] - 0.5)*psi_i[k] * sin(param.omega*(param.tau+param.dTau));

		for (k = 0; k < param.N+1; k++)
			psi_r[k] = psi_r[k] + H_i[k] * param.dTau / 2;

		for (k = 1; k < param.N; k++)
			H_r[k] = -(psi_r[k + 1] + psi_r[k - 1] - 2.0*psi_r[k]) / (2.0*pow(param.dx, 2))
			+ param.kappa*(x[k] - 0.5)*psi_r[k] * sin(param.omega*(param.tau + param.dTau));

		licznik += param.dTau;
	}
	fclose(rhoF);
	fclose(outF);

	return 0;
}