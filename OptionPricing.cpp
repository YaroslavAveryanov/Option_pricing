
#include "OptionPricing.h"
#include "Random.cpp"
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>

//Function definition
VanillaCall::VanillaCall(double& Spot_, double& Strike_, double& Rate_, double& Vol_, double& Time_, double& q_)
{
	Rate = Rate_ / 100.0;
	Vol = Vol_ / 100.0;
	q = q_ / 100.0;
	Spot = Spot_;
	Time = Time_;
	Strike = Strike_;
}

double VanillaCall::getPrice() const
{
	double d1 = 0.0;
	double d2 = 0.0;
	double price = 0.0;
	d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	d2 = d1 - Vol * sqrt(Time);
	price = Spot * CDF(d1) * exp(-q * Time) - Strike * exp(-Rate * Time) * CDF(d2);
	
	return price;
}

double VanillaCall::getDelta() const
{
	double d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	double result = CDF(d1) * exp(-q * Time);
	
	return result;
}

double VanillaCall::getGamma() const
{
	double d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	double result = PDF(d1) / (Spot * Vol * sqrt(Time)) * exp(-q * Time);
	
	return result;
}

/*
double VanillaCall::getVega() const
{
	double d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	double result = exp(-q * Time) * Spot * PDF(d1) * sqrt(Time);
	
	return result;
}
*/

double VanillaCall::getTheta() const
{
	double d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	double d2 = d1 - Vol * sqrt(Time);
	double result = -exp(-q * Time) * Spot * PDF(d1) * Vol / (2.0 * sqrt(Time)) - Rate * Strike * exp(-Rate * Time) * CDF(d2) + q * Spot * exp(-q * Time) * CDF(d1);
	
	return result;
}

VanillaPut::VanillaPut(double& Spot_, double& Strike_, double& Rate_, double& Vol_, double& Time_, double& q_)
{
	Rate = Rate_ / 100.0;
	Vol = Vol_ / 100.0;
	q = q_ / 100.0;
	Spot = Spot_;
	Time = Time_;
	Strike = Strike_;
}

double VanillaPut::getPrice() const
{
	double d1 = 0.0;
	double d2 = 0.0;
	double price = 0.0;
	
	d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	d2 = d1 - Vol * sqrt(Time);
	price = exp(-Rate * Time) * Strike * CDF(-d2) - Spot * exp(-q * Time) * CDF(-d1);
	
	return price;
}


double VanillaPut::getDelta() const
{
	double d1 = 0.0;
	double result;
	
	d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	
	result = exp(-q * Time) * (CDF(d1) - 1);
	
	return result;
}

/*
double VanillaPut::getVega() const
{
	double d1 = 0.0;
	double result;
	
	d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	result = Spot * exp(-q * Time) * PDF(d1) * sqrt(Time);
	
	return result;
}
*/

double VanillaPut::getTheta() const
{
	double d1 = 0.0;
	double d2 = 0.0;
	double result;
	
	d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	d2 = d1 - Vol * sqrt(Time);
	result = -exp(-q * Time) * Spot * PDF(d1) * Vol / (2.0 * sqrt(Time)) + Rate * Strike * exp(-Rate * Time) * CDF(-d2) - q * Spot * exp(-q * Time) * CDF(-d1);
	
	return result;
}


double VanillaPut::getGamma() const
{
	double d1 = 0.0;
	double result;
	
	d1 = (log(Spot / Strike) + (Rate - q + Vol * Vol / 2.0) * Time) / (Vol * sqrt(Time));
	result = exp(-q * Time) * PDF(d1) / (Spot * Vol * sqrt(Time));
	
	return result;
}

/*
AmericanCall::AmericanCall(double& Spot_, double& Strike_, double& Rate_, double& Vol_, double& Time_)
{
	Rate = Rate_ / 100.0;
	Vol = Vol_ / 100.0;
	Spot = Spot_;
	Time = Time_;
	Strike = Strike_;
	
	SpotPosition = 1;
	TimeSteps = 1000;
	deltaT = deltaX = 1;
	UnderlyingSteps = 1000;
}

double AmericanCall::getPrice() const 
{	
	//Can we build a Grid inside this function?
	//Find corresponding price using linear interpolation
	double Price;
	Price = Grid[0][SpotPosition - 1] + (Grid[0][SpotPosition] - Grid[0][SpotPosition - 1]) * (Spot - SpotArray[SpotPosition - 1]) / (SpotArray[SpotPosition] - SpotArray[SpotPosition - 1]);
	
	return Price;
}

double AmericanCall::getDelta() const 
{
	double resultLeft, resultRight, result;
	
	//How to calculate derivative with respect to S?
	resultLeft = (Grid[1][(SpotPosition - 1) + 1] - Grid[1][(SpotPosition - 1) - 1]) / 2.0 / deltaX / Spot;
	resultRight = (Grid[1][SpotPosition + 1] - Grid[1][SpotPosition - 1]) / 2.0 / deltaX / Spot;
	result = resultLeft + (resultRight - resultLeft) * (Spot - SpotArray[SpotPosition - 1]) / (SpotArray[SpotPosition] - SpotArray[SpotPosition - 1]);
	
	return result;
}

double AmericanCall::getGamma() const
{
	double resultLeft, resultRight, result;
	
	//What about the second derivative?
	resultLeft = (Grid[1][(SpotPosition - 1) + 1] - 2 * Grid[1][(SpotPosition - 1)] + Grid[1][(SpotPosition - 1) - 1]) / deltaX / deltaX / Spot / Spot - (Grid[1][(SpotPosition - 1) + 1] - Grid[1][(SpotPosition - 1) - 1]) / 2.0 / deltaX / Spot / Spot;
	resultRight = (Grid[1][SpotPosition + 1] - 2 * Grid[1][SpotPosition] + Grid[1][SpotPosition - 1]) / deltaX / deltaX / Spot / Spot - (Grid[1][SpotPosition + 1] - Grid[1][SpotPosition - 1]) / 2.0 / deltaX / Spot / Spot;
	result = resultLeft + (resultRight - resultLeft) * (Spot - SpotArray[SpotPosition - 1]) / (SpotArray[SpotPosition] - SpotArray[SpotPosition - 1]);	
	
	return result;
}

//Alternative solution? Disadvantages of the following code
double AmericanCall::getVega() const
{
	double deltaVol = 1;
	double new_Rate = Rate * 100.0;
	double new_Vol = Vol * 100.0 + deltaVol;
	AmericanCall* dummy = new AmericanCall((double&) Spot, (double&) Strike, new_Rate, new_Vol, (double&) Time);
	double result = (dummy->getPrice() - getPrice()) / (deltaVol / 100.0);
	
	delete dummy;
	return result;
}

double AmericanCall::getTheta() const
{
	double deltaVol;
	deltaVol = 0.15;
	
	return deltaVol;
}
*/

AmericanDifference::AmericanDifference(double& Spot_, double& Strike_, double& Rate_, double& Vol_, double& Time_, double& q_, std::string& Type_)
{
	Rate = Rate_ / 100.0;
	Vol = Vol_ / 100.0;
	Spot = Spot_;
	Time = Time_;
	Strike = Strike_;
	q = q_ / 100.0;
	Type = Type_;
	GridSize = 500;
	
	
	double dS = Spot / GridSize;
	int sGridSize = (int) (Strike / dS) * 2;
	double dt = pow(dS, 2) / (pow(Vol, 2) * pow(Strike, 2) * 4);
	int tGridSize = (int) (Time / dt) + 1;
	dt = Time / tGridSize;
	double interest = 1 / (1 + Rate * dt);
	
	int flag = Type == "Call" ? 1 : -1;
	
	spots = new double[sGridSize + 1];
	
	//double[][] values = new double[tGridSize + 1][sGridSize + 2];
	values = new double*[tGridSize + 1];
	
	for(int i = 0; i <= tGridSize; i++)
	{
		values[i] = new double[sGridSize + 2];
	}
	
	
	for(int i = 0; i <= sGridSize; i++)
	{
		spots[i] = dS * i;
		values[tGridSize][i] = std::max((double) 0, flag * (spots[i] - Strike));
	}
	
	for(int j = tGridSize - 1; j >= 0; j--)
	{
		for(int i = 1; i < sGridSize; i++)
		{
			double pu = 0.5 * (pow(Vol * i, 2) + (Rate - q) * i) * dt;
			double pd = 0.5 * (pow(Vol * i, 2) - (Rate - q) * i) * dt;
			double pm = 1 - pu - pd;
			values[j][i] = interest * (pu * values[j + 1][i + 1] + pm * values[j + 1][i] + pd * values[j + 1][i - 1]);
			values[j][i] = std::max(flag * (spots[i] - Strike), values[j][i]);
			values[j][0] = flag == 1 ? spots[sGridSize] - Strike : 0;
			values[j][sGridSize] = flag == 1 ? spots[sGridSize] - Strike : 0;
		}
	}
	
	output = new double[4];
	int sGridPoint = (int) (Strike / dS);
	output[0] = values[0][sGridPoint]; //price
    output[1] = (values[0][sGridPoint + 1] - values[0][sGridPoint - 1]) / (2 * dS); //delta
    output[2] = (values[0][sGridPoint + 1] - 2 * values[0][sGridPoint] + values[0][sGridPoint - 1])/pow(dS, 2); //gamma
    output[3] = (values[1][sGridPoint] - values[0][sGridPoint])/dt; //theta
}

AmericanDifference::~AmericanDifference()
{
	delete[] values;
	delete spots;
	delete output;
}

double AmericanDifference::getPrice() const 
{	
	return output[0];
}

double AmericanDifference::getDelta() const 
{
	return output[1];
}

double AmericanDifference::getGamma() const
{
	return output[2];
}

double AmericanDifference::getTheta() const
{
	return output[3];
}

AmericanBinomial::AmericanBinomial(double& Spot_, double& Strike_, double& Rate_, double& Vol_, double& Time_, double& q_, std::string& Type_)
{
	Rate = Rate_ / 100.0;
	Vol = Vol_ / 100.0;
	Spot = Spot_;
	Time = Time_;
	Strike = Strike_;
	q = q_ / 100.0;
	Type = Type_;
	TreeDepth = 500;
	
	double dt = Time / TreeDepth;
	double u = exp(Vol * sqrt(dt));
	double d = 1 / u;
	double p = (exp((Rate - q) * dt) - d) / (u - d);
	double interest = exp(-Rate * dt);
	double temp = 0;
	
	output = new double[4];
	int flag = Type == "Call" ? 1 : -1;
	values = new double[TreeDepth + 1];
	
	for(int i = 0; i <= TreeDepth; i++)
	{
		values[i] = std::max((double) 0, flag * (Spot * pow(u, i) * pow(d, TreeDepth - i) - Strike));
	}
	
	for(int j = TreeDepth - 1; j >= 0; j--)
	{
		for(int i = 0; i <= j; i++)
		{
			values[i] = interest * (p * values[i + 1] + (1 - p) * values[i]);
			values[i] = std::max(flag * (Spot * pow(u, i) * pow(d, j - i) - Strike), values[i]);
		}
		if (j == 2)
		{
			output[2] = ((values[2] - values[1]) / (Spot * pow(u, 2) - Spot) - (values[1] - values[0]) / (Spot - Spot * pow(d, 2))) / (0.5 * Spot * (pow(u, 2) - pow(d, 2)));
			temp = values[1];
		}
		if (j == 1)
		{
			output[1] = (values[1] - values[0]) / (Spot * (u - d));
		}
	}
	output[0] = values[0];
	output[3] = (temp - values[0]) / (2 * dt);
}	

AmericanBinomial::~AmericanBinomial()
{
	delete values;
	delete output;
}

double AmericanBinomial::getPrice() const 
{	
	return output[0];
}

double AmericanBinomial::getDelta() const 
{
	return output[1];
}

double AmericanBinomial::getGamma() const
{
	return output[2];
}

double AmericanBinomial::getTheta() const
{
	return output[3];
}

AmericanTrinomial::AmericanTrinomial(double& Spot_, double& Strike_, double& Rate_, double& Vol_, double& Time_, double& q_, std::string& Type_)
{
	Rate = Rate_ / 100.0;
	Vol = Vol_ / 100.0;
	Spot = Spot_;
	Time = Time_;
	Strike = Strike_;
	q = q_ / 100.0;
	Type = Type_;
	TreeDepth = 500;
	
	double dt = Time / TreeDepth;
	double u = exp(Vol * sqrt(2 * dt));
	double d = 1 / u;

	double pu = pow((exp((Rate - q) * dt / 2) - exp(-Vol * sqrt(dt / 2))) / (exp(Vol * sqrt(dt / 2)) - exp(-Vol * sqrt(dt / 2))), 2);
	double pd = pow((exp(Vol * sqrt(dt / 2)) - exp((Rate - q) * dt / 2)) / (exp(Vol * sqrt(dt / 2)) - exp(-Vol * sqrt(dt / 2))), 2);
	double pm = 1 - pu - pd;
	
	double interest = exp(-Rate * dt);
	output = new double[4];
	
	int flag = Type == "Call" ? 1 : -1;
	
	values = new double[2 * TreeDepth + 1];
	for(int i = 0; i <= 2 * TreeDepth; i++)
	{
		values[i] = std::max((double) 0, flag * (Spot * pow(u, std::max(i - TreeDepth, 0)) * pow(d, std::max(TreeDepth - i, 0)) - Strike)); 
	}
		
	for(int j = TreeDepth - 1; j >= 0; j--)
	{
		for(int i = 0; i <= 2 * j; i++)
		{
			values[i] = interest * (pu * values[i + 2] + pm * values[i + 1] + pd * values[i]); 
			values[i] = std::max(flag * (Spot * pow(u, std::max(i - j, 0)) * pow(d, std::max(j - i, 0)) - Strike), values[i]);
		}
		if (j == 1)
		{
			output[2] = ((values[2] - values[1]) / (Spot * u - Spot) - (values[1] - values[0]) / (Spot - Spot * d)) / (0.5 * Spot * (u - d));
			output[3] = values[1];
			output[1] = (values[2] - values[0]) / (Spot * (u - d));
		}
		
	}
	output[0] = values[0];
	output[3] = (output[3] - values[0]) / dt;
}

AmericanTrinomial::~AmericanTrinomial()
{
	delete values;
	delete output;
}

double AmericanTrinomial::getPrice() const 
{	
	return output[0];
}

double AmericanTrinomial::getDelta() const 
{
	return output[1];
}

double AmericanTrinomial::getGamma() const
{
	return output[2];
}

double AmericanTrinomial::getTheta() const
{
	return output[3];
}

