#include <iostream>
#include <cmath>
#include<boost/math/distributions.hpp>

using namespace std;

double erf(double x)
{  
	double y = 1.0 / ( 1.0 + 0.3275911 * x);   
	return 1 - (((((+ 1.061405429  * y - 1.453152027) * y + 1.421413741) * y - 0.284496736) * y + 0.254829592) * y) * exp(-x * x);      
}

double pdf(double x, double mu, double std_dev)
{
	static const double pi = 3.14159265; 
	return exp( -1 * (x - mu) * (x - mu) / (2 * std_dev * std_dev)) / (std_dev * sqrt(2 * pi));
}

double cdf(double x, double mu, double std_dev)
{
	return 0.5 * (1 + erf((x - mu) / (std_dev * sqrt(2.))));
}

std::pair<double,double> t_distribution_confidence_interval(const std::vector<double>& data, double alpha);

std::vector<std::vector<double>> kloops(std::vector<std::vector<double>> d, double alpha, int d_size){
    for (int i = 0; i<d_size; i++){
        for (int j = 0; j<d_size; j++){
            d[i][j] = alpha/(alpha + d[i][j]);
        }
    }
    return d;
}

	void UnivariateDistribution<T>::QuantileFunction(const std::vector<double> &p, std::vector<double> &y)
    int size = std::min(p.size(), y.size());
    for (int i = 0; i != size; ++i)
        y[i] = QuantileFunction(p[i]);
}
	boost::math::normal norm;
std::cout << QuantileFunction(complement(norm, 0.05)) << std::endl;	
}

int main()
{	
	double x, mu, std_dev;
	
	cout << "Insira os valores, na sequencia, de: x, mu, std_dev: ";
	cin >> x >> mu >> std_dev;

	cout << "PDF de x e: " << pdf(x, mu, std_dev) << endl;
	cout << "CDF de x e: " << cdf(x, mu, std_dev) << endl;
	cout << "CDF of x e (mu = 0, std_dev = 1): " << cdf(x, 0, 1) << endl;

	return 0;
}
