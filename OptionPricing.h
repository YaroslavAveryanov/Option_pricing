#ifndef OPTIONPRICING_H
#define OPTIONPRICING_H

#include <vector>
#include <string>

//Abstract class for Option - defined interface, overload payoff function in inherited classes
class Option
{
public:
	Option(){};
	virtual double getPrice() const = 0;
	virtual double getDelta() const = 0;
	virtual double getGamma() const = 0;
	//virtual double getVega() const = 0;
	virtual double getTheta() const = 0;
};

class VanillaCall : public Option
{
public:
	VanillaCall(double&, double&, double&, double&, double&, double&);
	virtual double getPrice() const;
	virtual double getDelta() const;
	virtual double getGamma() const;
	//virtual double getVega() const;
	virtual double getTheta() const;
private:
	double Spot;
	double Strike;
	double q; //annual dividend yield
	double Rate; //in % annualized
	double Vol; //in % annualized
	double Time; //time to maturity in years
};

class VanillaPut : public Option
{
public:
	VanillaPut(double&, double&, double&, double&, double&, double&);
	virtual double getPrice() const;
	virtual double getDelta() const;
	virtual double getGamma() const;
	//virtual double getVega() const;
	virtual double getTheta() const;
private:
	double Spot;
	double Strike;
	double q;
	double Rate; 
	double Vol; 
	double Time; 
};

class American : public Option
{
public:
		std::string Type;
		double* output;
};


class AmericanDifference : public American
{	
public:
	AmericanDifference(double&, double&, double&, double&, double&, double&, std::string&);
	virtual double getPrice() const;
	virtual double getDelta() const;
	virtual double getGamma() const;
	//virtual double getVega() const;
	virtual double getTheta() const;
	
	~AmericanDifference();
private:
	double Spot;
	double Strike;
	double Rate; //in % annualized
	double Vol; //in % annualized
	double Time; //time to maturity in years
	double q;
	
	int GridSize;
	
	double** values; //pricing grid
	double* spots; //possible spot values
};

class AmericanBinomial : public American
{
public:
	AmericanBinomial(double&, double&, double&, double&, double&, double&, std::string&);
	virtual double getPrice() const;
	virtual double getDelta() const;
	virtual double getGamma() const;
	//virtual double getVega() const;
	virtual double getTheta() const;
	
	~AmericanBinomial();
private:
	double Spot;
	double Strike;
	double Rate; //in % annualized
	double Vol; //in % annualized
	double Time; //time to maturity in years
	double q;
	
	int TreeDepth;
	double* values;
};

class AmericanTrinomial : public American
{
public:
	AmericanTrinomial(double&, double&, double&, double&, double&, double&, std::string&);
	virtual double getPrice() const;
	virtual double getDelta() const;
	virtual double getGamma() const;
	//virtual double getVega() const;
	virtual double getTheta() const;
	
	~AmericanTrinomial();
private:
	double Spot;
	double Strike;
	double Rate; //in % annualized
	double Vol; //in % annualized
	double Time; //time to maturity in years
	double q;
	
	int TreeDepth;
	double* values;
};

#endif
