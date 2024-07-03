#ifndef _MyModel_
#define _MyModel_

#include <vector>
#include "DNest4/code/DNest4.h"
#include "MyDistribution.h"
#include "Eigen/Dense"
#include "Eigen/Cholesky"

class MyModel
{
	private:
		DNest4::RJObject<MyDistribution> objects;

		double background;
		double extra_sigma; // Noise scale parameter

		// Parameters for the quasi-periodic extra noise
		double eta1, eta2, eta3, eta4;

		// The signal
		std::vector<long double> mu;
		void calculate_mu();

		// The covariance matrix for the data
		Eigen::MatrixXd C;
		void calculate_C();

		unsigned int staleness;

	public:
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

