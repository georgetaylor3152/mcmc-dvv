"""
A class that defines a Markov Chain Monte Carlo
inversion
"""

import numpy as np

from scipy.sparse import csr_matrix
from scipy.stats import norm

class MarkovChainMonteCarlo:
    def __init__(self, observations, g, weights, mlen, prior_low,
                 prior_high):
        self.observations = observations
        self.g = g
        self.weights = weights # Data weights (1 / variance)
        self.mlen = mlen # Number of correlation functions
        # The prior bounds for uniform prior
        self.prior_low = prior_low
        self.prior_high = prior_high

    def generate_prior(self):
        """
        Generate a prior starting model from a uniform probability
        distribution
        """
        prior_model = np.random.uniform(low=self.prior_low, 
                                        high=self.prior_high, 
                                        size=self.mlen)

        return prior_model

    def gaussian_mtrial(self, mcurr, std):
        """
        Generates a proposal model from a small perturbation
        of the current model, mcurr. The perturbation is drawn from
        a Gaussian distribution with zero mean and standard
        deviation = std
        """
        # Generate the model perturbation
        perturb = np.random.normal(0.0, std, len(mcurr))
        # Perturb the current model
        mtrial = mcurr + perturb

        return mtrial

    def predict_observations(self, m):
        """
        Calculates the set of observations that are predicted for model m,
        given the physics of the problem in g. d_predict = G * m
        """
        d_predicted = self.g.dot(m)
        return d_predicted

    def evaluate_likelihood(self, mtrial):
        """
        Evaluates the likelihood of our observations, given a
        proposed model, mtrial
        """
        # Find the predicted observations for mtrial
        d_predicted = self.predict_observations(mtrial)
        # Calculate the vector of residuals between the predicted and
        # actual observations
        residuals = self.observations - d_predicted
        # Make residuals a "sparse" matrix so we can multiply with cd_inv
        exponent = np.dot((residuals * self.weights),residuals)
        # Calculate the likelihood
        likelihood = -0.5 * exponent

        return likelihood

    def calculate_prior_uniform(self, mtrial):
        """
        Calculates whether the prposed model violates the uniform
        prior. Returns 1 when it does not, 0 when it does.
        """
        if np.any([mtrial > self.prior_high]) or np.any([mtrial < self.prior_low]):
            probability = 0.
        else:
            probability = 1.
        return probability

    def sensitivity_test(self, accept_reject, proposal_std): 
        """
        Performs a sensitivity test to determine a new step size for the 
        proposals
        """
        # Rejection ratio
        no_rejects = len(accept_reject[np.where(accept_reject == False)])
        rejection_ratio = no_rejects / 100.
        print("The rejection ratio since the last test is: {}".format(rejection_ratio))

        if rejection_ratio > 0.73 and rejection_ratio < 0.79:
            new_proposal = proposal_std
        else:
            # We want to accept 23.4% of proposals
            difference = 0.766 / rejection_ratio
            new_proposal = difference * proposal_std

        print("The new step size is: {}".format(new_proposal))

        return new_proposal

    def do_mcmc(self, iterations, proposal_std):
        """ 
        Performs the MCMC sampling of the posterior probability
        """
        # Generate a model from the prior as an initial guess
        mcurr = self.generate_prior()
        # Remove the mean of the model
        mcurr -= np.mean(mcurr)

        # Array for storing the posterior distribution, posterior likelihoods,
        # and whether the trial was accepted/rejected
        posterior_distribution = np.zeros((self.mlen,iterations))
        posterior_likelihoods = np.zeros((iterations))
        accept_reject = np.zeros((100))
        for i in range(0, iterations):
            # Calculate the likelihood of the current model satisfying the data
            current_likelihood = self.evaluate_likelihood(mcurr)
            # Calculate the probability of the current model
            current_prob = current_likelihood

            # Calculate a new trial model
            mtrial = self.gaussian_mtrial(mcurr, proposal_std)
            # Remove the mean of the model
            mtrial -= np.mean(mtrial)
            # Check the prior probability of the proposed model
            prior_prob_trial = self.calculate_prior_uniform(mtrial)

            if not prior_prob_trial:
                # The trial is rejected
                accept_reject[i % 100] = False
                posterior_likelihoods[i] = current_likelihood
                posterior_distribution[:,i] = mcurr

                # Check for sensitivity test
                if not (i == 0) and (i % 100 == 0):
                    print(50*"-")
                    print("Sensitivity test for iteration number {}:".format(i))
                    proposal_std  = self.sensitivity_test(accept_reject, proposal_std)
                    print("The current likelihood is {}".format(
                        posterior_likelihoods[i]))
                    accept_reject = np.zeros((100))
                    print(50*"-")

                continue

            # Evaluate the likelihood of the trial model
            trial_likelihood = self.evaluate_likelihood(mtrial)
            # Calculate the probability of the trial model
            trial_prob = trial_likelihood
            
            # If the trial model has better or equal likelihood
            # than the current model, it becomes in the new
            # current model
            if trial_prob >= current_prob:
                #print("The likelihood is higher!")
                # The trial is accepted
                mcurr = mtrial
                accept_reject[i % 100] = True
                posterior_likelihoods[i] = trial_likelihood

            # Even if the likelihood is lower, we still may accept the
            # trial model. The probability that we do so is determined
            # the ratio of the two likelihoods
            else:
                # Likelihood ratio
                ratio = np.exp(trial_prob - current_prob)
                # Determine a random number between 0 and 1
                acceptance_thresh = np.random.rand(1)
                # If the likelihood ratio exceeds the acceptance threshold
                # We accept the new model
                if ratio > acceptance_thresh:
                    # The trial is accepted
                    mcurr = mtrial
                    accept_reject[i % 100] = True
                    posterior_likelihoods[i] = trial_likelihood
                else:
                    # The trial is rejected
                    #print("Proposal rejected due to likelihood")
                    accept_reject[i % 100] = False
                    posterior_likelihoods[i] = current_likelihood

            # Store the new model, whether the trial was accepted or rejected
            posterior_distribution[:,i] = mcurr

            if not (i == 0) and (i % 100 == 0):
                print(50*"-")
                print("Sensitivity test for iteration number {}:".format(i))
                proposal_std  = self.sensitivity_test(accept_reject, proposal_std)
                print("The current likelihood is {}".format(
                    posterior_likelihoods[i]))
                accept_reject = np.zeros((100))
                print(50*"-")

        return posterior_distribution, posterior_likelihoods

    def calculate_acception_rate(self, accept_reject):
        """
        Calculates the rate at which the MCMC run accepted/rejected models
        """
        # Store the acceptance ratio at each iteration
        acceptance_ratios = np.zeros(len(accept_reject))
        for i, iteration in enumerate(accept_reject):
            # Grab all results up to this iteration
            results = accept_reject[:i]
            # Calculate how many accepts up to this point
            no_accepts = len(results[np.where(results == True)])
            no_rejects = len(results[np.where(results == False)])
            # Must avoid devision by zero
            if no_rejects == 0:
                acceptance_ratios[i] = 1.
            else:
                # Ratio
                acceptance_ratios[i] = (len(results) - no_rejects) / len(results)

        return acceptance_ratios

        


