#include "MyDistribution.h"
#include <cmath>

MyDistribution::MyDistribution()
{

}

double MyDistribution::cauchy(double x, double center, double scaler)
{
	return center + scaler * tan(M_PI*(0.97*x - 0.485));
}

/// @brief Generate random sample from a cauchy distribution
/// @param rng Random number generator
/// @param center Center of the distribution
/// @param scaler Scaler of the distribution
/// @return sample
double MyDistribution::cauchy_sample(DNest4::RNG& rng, double center = 0, double scaler = 1) 
{
	return cauchy(rng.rand(), center, scaler);
}

/// @brief Generate random sample from a gaussian distribution
/// @param rng Random number generator
/// @param mean Expected value of the gaussian
/// @param sigma Std dev of the gaussian
/// @return sample
double MyDistribution::gaussian_sample(DNest4::RNG& rng, double mean = 0, double sigma = 1) 
{
	// generate gaussian sample with Box-Muller transform
	double u1 = rng.rand();
	double u2 = rng.rand();
	double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
	return z * sigma + mean;
}

/// @brief Randomly perturb a sample from a cauchy distribution, preserving distribution characteristics
/// @param rng Random number generator
/// @param sample Original sample
/// @param center Original center
/// @param scaler Original scaler
/// @return Perturbed sample
double MyDistribution::perturb_cauchy(DNest4::RNG& rng, double sample, double center = 0, double scaler = 1)
{
	// inverse cauchy to retrieve original random number
	double orig_rn = (atan((sample - center)/scaler)/M_PI + 0.485)/0.97;
	double new_rn = orig_rn + rng.randh();
	DNest4::wrap(new_rn, 0., 1.);
	return cauchy(new_rn, center, scaler);
}

/// @brief Randomly perturb a sample from a normal distribution, preserving distribution characteristics
/// @param rng Random number generator
/// @param sample Original sample
/// @param mean Original expected value
/// @param sigma Original std dev
/// @param perturb_sigma Additional std dev for the perturbation
/// @return Perturbed sample
double MyDistribution::perturb_gaussian(DNest4::RNG& rng, double sample, double mean = 0, double sigma = 1, double perturb_sigma = 1)
{
	double perturbation = gaussian_sample(rng, 0, perturb_sigma);
	double perturbed_sample = sample + perturbation;
	// ensure the same mu and sigma for the returned sample
	return mean + (perturbed_sample - mean) * (sigma / std::sqrt(pow(sigma, 2) + pow(perturb_sigma, 2)));
}

void MyDistribution::from_prior(DNest4::RNG& rng)
{
	// Cauchy prior centered on 5.901 = log(365 days).
	// center = 5.901 + tan(M_PI*(0.97*rng.rand() - 0.485));
	center = cauchy_sample(rng, 5.901);
	width = 0.1 + 2.9*rng.rand();
	mu = exp(cauchy_sample(rng));
}

double MyDistribution::perturb_hyperparameters(DNest4::RNG& rng)
{
	double logH = 0.;

	int which = rng.rand_int(3);

	if(which == 0)
	{
		center = perturb_cauchy(rng, center, 5.901);
	}
	else if(which == 1)
	{
		width += 2.9*rng.randh();
		DNest4::wrap(width, 0.1, 3.);
	}
	else
	{
		mu = log(mu);
		mu = perturb_cauchy(rng, mu);
		mu = exp(mu);
	}
	return logH;
}

// vec[0] = "position" (log-period)
// vec[1] = amplitude
// vec[2] = phase
// vec[3] = v0
// vec[4] = viewing angle

double MyDistribution::log_pdf(const std::vector<double>& vec) const
{
	if(vec[1] < 0. ||
			vec[2] < 0. || vec[2] > 2.*M_PI ||
			vec[3] < 0. || vec[3] > 0.8189776 ||
			vec[4] < 0. || vec[4] > 2.*M_PI)
		return -1E300;

	return  -log(2.*width) - abs(vec[0] - center)/width
		-log(mu) - vec[1]/mu
		+ 2.1*log(1. - vec[3]/0.995);
}

void MyDistribution::from_uniform(std::vector<double>& vec) const
{
	if(vec[0] < 0.5)
		vec[0] = center + width*log(2.*vec[0]);
	else
		vec[0] = center - width*log(2. - 2.*vec[0]);
	vec[1] = -mu*log(1. - vec[1]);
	vec[2] = 2.*M_PI*vec[2];
	vec[3] = 1. - pow(1. - 0.995*vec[3], 1./3.1);
	vec[4] = 2.*M_PI*vec[4];
}

void MyDistribution::to_uniform(std::vector<double>& vec) const
{
	if(vec[0] < center)
		vec[0] = 0.5*exp((vec[0] - center)/width);
	else
		vec[0] = 1. - 0.5*exp((center - vec[0])/width);
	vec[1] = 1. - exp(-vec[1]/mu);
	vec[2] = vec[2]/(2.*M_PI);
	vec[3] = 1. - pow(1. - vec[3]/0.995, 3.1);
	vec[4] = vec[4]/(2.*M_PI);
}

void MyDistribution::print(std::ostream& out) const
{
	out<<center<<' '<<width<<' '<<mu<<' ';
}

