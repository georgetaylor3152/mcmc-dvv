============================================================
Markov chain Monte Carlo code for calculating dv time series
============================================================

Usage
=====
The module is imported using:

``from mcmc_inversion import MarkovChainMonteCarlo``

An MCMC object must then be instantiated in your code:

``mcmc_inversion = MarkovChainMonteCarlo(data, g, weights, n_corrs, prior_low, prior_high)``

Where:
data = A numpy vector of your observations, "d" vector as described by Brenguier et al. (2014)
g = A 2D numpy array, the "G"  matrix of your problem, as described by Brenguier et al. (2014)
weights = A numpy vector of data weights, equal to 1/variances of your data vector.
n_corrs = An integer, the number of correlation functions you used, the dv time series
          will have the same number of points.
prior_low = A float, the lower bound of your prior expectation for the dv values
prior_high = A float, the upper bound of your prior expectation for the dv values

To begin the MCMC chain, run the following:

``distribution, likelihood = mcmc_inversion.do_mcmc(250000, 3.0e-6)``

The first argument is an integer, and is the number of iterations to run the chain for.
This value **must** be divisible by 100.

The second argument is a float, and is the initial size of the model perturbation. This
value will be continuously updated by the MCMC code every 100 iterations.

The MCMC algorithm returns three outputs:

The first is a 2D array of the posterior distribution. Your dv models are stored in the
**columns** of this array.

A vector containing the likelihood of each model in the posterior distribution.
length(likelihood) = number of columns in distribution

License
=======
This project is licensed under the MIT License - see the LICENSE file for details
