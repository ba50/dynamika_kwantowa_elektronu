#define PI 3.1416

#include <string>

struct Parameters{
	std::string IN_filename;
	std::string OUT_filename;
	std::string RHO_filename;

	unsigned int N;
	unsigned int n;

	double dx;
	double dTau;
	double kappa;
	double omega;
	double tau;
	double end;
};