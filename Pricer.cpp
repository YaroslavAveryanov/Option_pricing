
#include "OptionPricing.cpp"
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

int main()
{
	double q, Rate, Strike, Spot, Vol, Time;
	string Style;
	string Type;
	
	
	cout << "Dividend Yield = "; cin  >> q;
	cout << "Interest Rate = "; cin >> Rate;
	cout << "Strike = "; cin >> Strike;
	cout << "Spot = "; cin >> Spot;
	cout << "Volatility = "; cin >> Vol;
	cout << "Time = "; cin >> Time;
	cout << "Style = "; cin >> Style;
	cout << "Type = "; cin >> Type;
	
	if (Style == "European")
	{
		if (Type == "Call")
		{
			VanillaCall* Option = new VanillaCall(Spot, Strike, Rate, Vol, Time, q);
			cout << "Price of European call option = " << Option->getPrice() << endl;
			cout << "Delta of European call option = " << Option->getDelta() << endl;
			cout << "Theta of European call option = " << Option->getTheta() << endl;
			//cout << "Vega of European call option = " << Option->getVega() << endl;
			cout << "Gamma of European call option = " << Option->getGamma() << endl;
			delete Option;
		}
		else if (Type == "Put")
		{
			VanillaPut* Option = new VanillaPut(Spot, Strike, Rate, Vol, Time, q);
			cout << "Price of European put option = " << Option->getPrice() << endl;
			cout << "Delta of European put option = " << Option->getDelta() << endl;
			cout << "Theta of European put option = " << Option->getTheta() << endl;
			//cout << "Vega of European put option = " << Option->getVega() << endl;
			cout << "Gamma of European put option = " << Option->getGamma() << endl;
			delete Option;
		}
	}
	else if (Style == "American")
	{
		string American_Pricing;
		cout << "Method for pricing = "; cin >> American_Pricing;
		if (American_Pricing == "Finite_difference")
		{
			AmericanDifference* Option = new AmericanDifference(Spot, Strike, Rate, Vol, Time, q, Type);
			cout << "Price of American option = " << Option->getPrice() << endl;
			cout << "Delta of American option = " << Option->getDelta() << endl;
			cout << "Theta of American option = " << Option->getTheta() << endl;
			cout << "Gamma of American option = " << Option->getGamma() << endl;
			delete Option;
		}
		else if (American_Pricing == "Binomial")
		{
			AmericanBinomial* Option = new AmericanBinomial(Spot, Strike, Rate, Vol, Time, q, Type);
			cout << "Price of American option = " << Option->getPrice() << endl;
			cout << "Delta of American option = " << Option->getDelta() << endl;
			cout << "Theta of American option = " << Option->getTheta() << endl;
			cout << "Gamma of American option = " << Option->getGamma() << endl;
			delete Option;
		}
		else if (American_Pricing == "Trinomial")
		{
			AmericanTrinomial* Option = new AmericanTrinomial(Spot, Strike, Rate, Vol, Time, q, Type);
			cout << "Price of American option = " << Option->getPrice() << endl;
			cout << "Delta of American option = " << Option->getDelta() << endl;
			cout << "Theta of American option = " << Option->getTheta() << endl;
			cout << "Gamma of American option = " << Option->getGamma() << endl;
			delete Option;
		}
		else cout << "No such style of method for pricing" << endl;
	}
	else cout << "No such style of option" << endl;
	
	return 0;
}
