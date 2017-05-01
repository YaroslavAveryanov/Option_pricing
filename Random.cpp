
#include "Random.h"
#include <iostream>
#include <cmath>

//Available in other code files?
const double PI = 3.141592653589793238463;

double CDF(double x)
{
	double L = 0.0;
	double K = 0.0;
	double dCDF = 0.0;
	const double a1 = 0.31938153;
	const double a2 = -0.356563782;
	const double a3 = 1.781477937;
	const double a4 = -1.821255978;
	const double a5 = 1.330274429;
	L = std::abs(x);
	K = 1.0 / (1.0 + 0.2316419 * L);
	
	dCDF = 1.0 - 1.0 / sqrt(2 * PI) *
	exp(-L * L / 2.0) * (a1 * K + a2 * K * K + a3 * pow(K, 3.0) +
	a4 * pow(K, 4.0) + a5 * pow(K, 5.0));
	
	if (x < 0)
	{
		return 1.0 - dCDF;
	}
	else
	{
		return dCDF;
	}
	
}

double newCDF(double x)
{
	// constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}

double PDF(double x)
{
	double dPDF = 0.0;
	
	dPDF = 1.0 / (sqrt(2 * PI)) * exp(- x * x / 2);
	
	return dPDF;
}



