# Copyright 2012 Oliver Serang
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

from utilities import *
import numpy
import math
import sys
import getopt
import cProfile
import itertools
from copy import deepcopy
from pprint import pprint

from numerics import *

# note: performance can be improved by making functions that don't use
# "self" static
class Base:
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        self.load(filename, **Kwargs)
        self.ploidy_range = ploidy_range
        self.sigma_range = sigma_range

    def load(self, filename, heights_or_areas):
        try:
            f = open(filename)
            individuals_to_data = {}
            all_lines = f.readlines()
            if len(all_lines) == 0:
                raise Exception("cannot load an empty file")
            for line in all_lines:
                line_lst = line.split()
                if len(line_lst) == 5:
                    individual, x_height, y_height, x_area, y_area = line_lst
                    if heights_or_areas == 'heights':
                        pair = numpy.array( [float(x_height), float(y_height)] )
                    elif heights_or_areas == 'areas':
                        pair = numpy.array( [float(x_area), float(y_area)] )
                    else:
                        raise Exception('load_scatterplot needs to specify height or area')
                elif len(line_lst) == 3:
                    if heights_or_areas == 'areas':
                        raise Exception('heights_or_areas does not affect 3 column data sets')
                    individual, x, y = line_lst
                    pair = numpy.array([float(x), float(y)])
                else:
                    raise Exception('cannot load file with number of columns != 3 or 5')

                if numpy.linalg.norm(pair) > 0.0:
                    # add this to the data set
                    if individual in individuals_to_data:
                        individuals_to_data[ individual ].append( pair )
                    else:
                        individuals_to_data[ individual ] = [ pair ]
            self.individuals_to_data = individuals_to_data
        except Exception as err:
            raise Exception('Error loading file: ' + str(err))

    def get_mean_for_each_individual(self):
        individuals_to_mean_data = {}
        for individual, data in self.individuals_to_data.items():
            mean = sum( data ) / float(len(data))
            individuals_to_mean_data[individual] = mean
        return individuals_to_mean_data

    def get_individuals_to_log_likelihoods(self, ploidy, sigma, **Kwargs):
        # takes the product of all data points for each genotype state
        # for each individual
        individuals_to_log_likelihoods = {}
        for individual in self.individuals_to_data:
            sum_log_likelihoods = None
            for pair in self.individuals_to_data[individual]:
                log_likelihoods = self.log_likelihoods_of_genotype_states(ploidy, sigma, pair)

                if sum_log_likelihoods == None:
                    sum_log_likelihoods = log_likelihoods
                else:
                    for g in log_likelihoods:
                        sum_log_likelihoods[g] += log_likelihoods[g]
            individuals_to_log_likelihoods[individual] = sum_log_likelihoods
        return individuals_to_log_likelihoods

    def greedy(self, individuals_to_log_likelihoods, allowed_genotypes):
        individuals_to_best_configurations = {}
        overall_sum_log_likelihoods = 0.0
        for individual, log_likelihoods in individuals_to_log_likelihoods.items():
            best_log_likelihood, best_genotype = max([ (log_like, g) for g, log_like in log_likelihoods.items() if g in allowed_genotypes ])
            individuals_to_best_configurations[ individual ] = best_genotype
            overall_sum_log_likelihoods += best_log_likelihood
        return overall_sum_log_likelihoods, individuals_to_best_configurations

    def log_likelihood_of_configuration(self, inds_to_genos, individuals_to_log_likelihoods):
        return sum([ individuals_to_log_likelihoods[ind][geno] for ind, geno in inds_to_genos.items() ])

    def log_likelihood_of_genotype(self, sigma, x_and_y_G, x_and_y):
        ## fixme: this should to include the fact that the observed x
        ## and y values can never be negative
        theoretical_x_and_y_normalized = numpy.array(x_and_y_G) / float(sum(x_and_y_G))
        x_and_y_normalized = x_and_y / float(sum(x_and_y))
        
        dst = numpy.linalg.norm(x_and_y_normalized - theoretical_x_and_y_normalized)
        return log_normal_pdf(dst, sigma )

    def log_likelihoods_of_genotype_states(self, ploidy, sigma, x_and_y):
        dist = {}
        normalization_constant = log(0)
        for x_G in xrange(0, ploidy+1):
            y_G = ploidy - x_G
            geno = (x_G, y_G)
            log_like = self.log_likelihood_of_genotype(sigma, geno, x_and_y)
            dist[ geno ] = log_like
        return dist

    def optimal(self):
        pass
    def optimal_greedy(self):
        pass

    def best_log_score_posterior_configuration_parameters_for_MAPs(self, log_joints_and_inds_to_genotypes_and_parameters):
        # over all parameters and results, for each set of parameters
        # to marginalize, take take the MAP. then marginalize over the
        # parameters to marginalize to get a posterior estimate.

        # build a map of marg vars to results; then the MAP of these
        # results can be taken for each marg var config
        params_to_MAPs = {}
        best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars = []
        for log_joint, config, params in log_joints_and_inds_to_genotypes_and_parameters:
            marg_vars_and_outcomes = tuple( [ (param, val) for param, val in params.items() if param not in self.vars_for_MAP ] )
            MAP_vars_and_outcomes = [ (var, params[var]) for var in self.vars_for_MAP ]

            if marg_vars_and_outcomes in params_to_MAPs:
                params_to_MAPs[marg_vars_and_outcomes].append( (log_joint, MAP_vars_and_outcomes, config) )
            else:
                params_to_MAPs[marg_vars_and_outcomes] = [ (log_joint, MAP_vars_and_outcomes, config) ]

        for marg_vars_and_outcomes, results in params_to_MAPs.items():
            log_joint, MAP_vars_and_outcomes, config = max(results)
            params_dict = dict(marg_vars_and_outcomes)
            params_dict.update( dict(MAP_vars_and_outcomes) )
            best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars.append( (log_joint, config, params_dict) )

        log_joint_max, inds_to_genotypes_max, parameters_max = max(best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars)

        best_posterior = math.exp( log_joint_max - log_sum([log_joint for log_joint, inds_to_genotypes, params in best_log_joints_and_inds_to_genotypes_and_parameters_for_marg_vars]) )
        return log_joint_max, best_posterior, inds_to_genotypes_max, parameters_max

    def naive_high_quality_genotypes(self, naive_posterior_threshold, individuals_to_genotypes, all_parameter_map, prob_geno_fname):
        individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(**all_parameter_map)
        f = open(str(prob_geno_fname), 'w+')
        high_qual_inds_to_genos = {}
        flag = 0
        for ind, log_likelihoods in individuals_to_log_likelihoods.items():
            flag = flag + 1
            log_total = log_sum( log_likelihoods.values() )
            log_chosen = log_likelihoods[ individuals_to_genotypes[ind] ]
            if math.exp(log_chosen - log_total) > naive_posterior_threshold:
                high_qual_inds_to_genos[ind] = individuals_to_genotypes[ind]
            if flag == 1:
               f.write(str(log_likelihoods.keys()))
               f.write('\n')
            f.write(ind)
            f.write( ' ' )
            s = str(log_likelihoods.values())
            f.write(s)
            f.write('\n')
        f.close()    
        return high_qual_inds_to_genos

class SharedPloidyInference(Base):
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        Base.__init__(self, filename, ploidy_range, sigma_range, **Kwargs)
        self.vars_for_MAP = ['sigma']

    def log_parameter_prior(self, ploidy):
        # treats genotypes as uniform given ploidy
        return log(1/float(ploidy+1))

    def optimal(self, display_progress):
        return self.optimal_greedy(display_progress)

    def get_inds_to_genotypes_given_parent_parameters(self, parameter, sigma, **Kwargs):
        geno_set = set(parameter)
        inds_to_log_likes = self.get_individuals_to_log_likelihoods(sigma=sigma, **Kwargs)
        # return the configuration, not the score
        return self.greedy(inds_to_log_likes, geno_set)[1]

    def optimal_greedy(self, display_progress):
        log_joints_and_inds_to_genotypes_and_parameters = []
        for sigma in self.sigma_range:
            for ploidy in self.ploidy_range:
                individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(ploidy, sigma)
                log_pr_D_given_G, inds_to_genotypes = self.greedy(individuals_to_log_likelihoods, set([ (i, ploidy-i) for i in xrange(0, ploidy+1) ]) )
                log_joint = log_pr_D_given_G + self.log_parameter_prior(ploidy)
                log_joints_and_inds_to_genotypes_and_parameters.append((log_pr_D_given_G, inds_to_genotypes,  {'sigma':sigma, 'ploidy':ploidy}))
                if display_progress:
                    print sigma, ploidy, log_joint
        return self.best_log_score_posterior_configuration_parameters_for_MAPs(log_joints_and_inds_to_genotypes_and_parameters)

# rewrite functors with two inherited classes from here
class PopulationInference(Base):
    def __init__(self, filename, ploidy_range, sigma_range, vars_for_MAP, **Kwargs):
        Base.__init__(self, filename, ploidy_range, sigma_range, **Kwargs)
        self.vars_for_MAP = set(vars_for_MAP)
    def optimal_greedy(self, display_progress):
        log_joints_and_inds_to_genotypes_and_parameters = []
        for sigma in self.sigma_range:
            for ploidy in self.ploidy_range:
                individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(ploidy, sigma)
                for parameter in self.parameter_range_generator(ploidy = ploidy):
                    T = self.T_generator(ploidy, parameter)
                    T_nonzero = nonzero_distribution(T)
                    log_joint, individuals_to_genotypes = self.single_greedy(ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, sigma)
                    log_joints_and_inds_to_genotypes_and_parameters.append((log_joint, individuals_to_genotypes, {'sigma':sigma, 'ploidy':ploidy, 'parameter':parameter}))
                    if display_progress:
                        print sigma, ploidy, parameter, log_joint
                        #log_total = log_likelihoods.values()
                        #print log_total

        return self.best_log_score_posterior_configuration_parameters_for_MAPs(log_joints_and_inds_to_genotypes_and_parameters)

    def single_greedy(self, ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, sigma):
        # get the greedy result for this set of parameters
        greedy_log_pr_D_given_G, greedy_inds_to_genotypes  = self.greedy(individuals_to_log_likelihoods, set(T_nonzero))
        greedy_C = histogram(greedy_inds_to_genotypes)
        greedy_log_pr_C_given_T = log_probability_of_histogram_given_distribution(greedy_C, T_nonzero)
        greedy_log_pr_D_and_C_given_G_and_T = greedy_log_pr_D_given_G + greedy_log_pr_C_given_T
        log_prior = self.log_parameter_prior(ploidy, parameter=parameter, sigma=sigma)
        greedy_log_score = greedy_log_pr_D_and_C_given_G_and_T + log_prior
        return greedy_log_score, greedy_inds_to_genotypes

    def optimal(self, display_progress, epsilon = 0.01):
        if display_progress:
            print 'seeding with greedy search'
        best_log_score, posterior, config, params = self.optimal_greedy(display_progress)
        if display_progress:
            print 'starting optimal search'

        if self.number_of_marginalized_configs > 0:
            threshold_buffer_for_marg = log(epsilon / float(self.number_of_marginalized_configs))
        else:
            # there is only one confuguration to try, so the buffer
            # can be 0
            threshold_buffer_for_marg = 0

        log_joints_and_inds_to_genotypes_and_parameters = []
        for sigma in self.sigma_range:
            for ploidy in self.ploidy_range:
                individuals_to_log_likelihoods = self.get_individuals_to_log_likelihoods(ploidy, sigma)
                for parameter in self.parameter_range_generator(ploidy = ploidy):
                    T = self.T_generator(ploidy, parameter)
                    T_nonzero = nonzero_distribution(T)

                    # note: this unecessarily sorts the individuals
                    # each iteration. it is arranged this way for
                    # simplicity; plus the performance gain would be
                    # small since runtime is dominated by branch and
                    # bound
                    ordered_individuals, ordered_genotypes = self.get_ordered_individuals_and_genotypes(T_nonzero)

                    # get the greedy log joint and configuration
                    greedy_log_joint, greedy_individuals_to_genotypes = self.single_greedy(ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, sigma)

                    # get a different greedy result based on the best
                    # histogram C
                    hist_greedy_log_joint, hist_greedy_individuals_to_genotypes = self.single_histogram_greedy(ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, ordered_individuals, ordered_genotypes, sigma)

                    # use the greedy results to seed the branch and
                    # bound. also, if a known configuration is
                    # sufficiently better, numerically bound.
                    # furthermore, give a small margin of 0.05 to
                    # prevent bounding from numerical error from
                    # preventing any solutions (in case the optimal
                    # isn't much better than the greedy)

                    # since the branch and bound doesn't use prior,
                    # and all of these do use a prior, then subtract
                    # the prior for this set of parameters from the
                    # threshold; it will be the same as adding it into
                    # the branch and bound. we want the set of
                    # parameters with maximum branch and bound result
                    # (including parameter prior)
                    log_prior = self.log_parameter_prior(ploidy, parameter=parameter, sigma=sigma)
                    threshold = max(greedy_log_joint, hist_greedy_log_joint, best_log_score + threshold_buffer_for_marg) - log_prior - 0.05

                    log_joint, individuals_to_genotypes = self.count_branch_and_bound(individuals_to_log_likelihoods, T_nonzero, threshold, ordered_individuals, ordered_genotypes)
                    log_joint += log_prior

                    best_log_score = max(best_log_score, log_joint)

                    if display_progress:
                        print sigma, ploidy, parameter, log_joint

                    log_joints_and_inds_to_genotypes_and_parameters.append((log_joint, individuals_to_genotypes, {'sigma':sigma, 'ploidy':ploidy, 'parameter':parameter}))
        return self.best_log_score_posterior_configuration_parameters_for_MAPs(log_joints_and_inds_to_genotypes_and_parameters)

    def single_histogram_greedy(self, ploidy, parameter, T_nonzero, individuals_to_log_likelihoods, ordered_individuals, ordered_genotypes, sigma):
        hist_greedy_C = self.greedy_best_possible_histogram_for_theoretical(T_nonzero, len(individuals_to_log_likelihoods))
        hist_greedy_inds_to_genotypes = self.get_individuals_to_genotypes_from_histogram(ordered_individuals, ordered_genotypes, hist_greedy_C)
        hist_greedy_log_pr_C_given_T = log_probability_of_histogram_given_distribution(hist_greedy_C, T_nonzero)
        hist_greedy_log_pr_D_given_G = sum([individuals_to_log_likelihoods[ind][geno] for ind, geno in hist_greedy_inds_to_genotypes.items() ])
        log_prior = self.log_parameter_prior(ploidy, parameter = parameter, sigma = sigma)
        hist_greedy_log_score = hist_greedy_log_pr_D_given_G + hist_greedy_log_pr_C_given_T + log_prior
        return hist_greedy_log_score, hist_greedy_inds_to_genotypes

    def greedy_best_possible_histogram_for_theoretical(self, T_nonzero, n):
        C = {}
        # use sorted order so that results aren't stochastic; I don't
        # believe dictionary items are not guaranteed to be in the
        # same order
        for i, g in enumerate(sorted(T_nonzero)):
            if i == len(T_nonzero)-1:
                break
            expected = int(n*T_nonzero[g])
            C[g] = expected
        # set the last one (it's the one that broke the loop)
        C[g] = n - sum(C.values())
        return C

    def get_ordered_individuals_and_genotypes(self, T_nonzero):
        # sort the genotypes lexicographically (after L1 normalization)
        ordered_individuals = [ (tuple(pair/sum(pair)), individual) for individual, pair in self.get_mean_for_each_individual().items() ]
        ordered_individuals = sorted(ordered_individuals)[::-1]
        # sort the theoretical genotypes lexicographically
        ordered_genotypes = sorted( list(T_nonzero) )[::-1]
        return ordered_individuals, ordered_genotypes

    def count_branch_and_bound(self, individuals_to_log_likelihoods, T_nonzero, threshold, ordered_individuals, ordered_genotypes):
        # arrange the log_likelihoods in the same order
        geometrically_sorted_individuals_and_log_likelihoods = [ (individual, individuals_to_log_likelihoods[individual]) for value, individual in ordered_individuals ]

        best_score, genotypes_to_counts = self.count_branch_and_bound_recursive(geometrically_sorted_individuals_and_log_likelihoods, T_nonzero, ordered_genotypes, {}, 0, 0.0, threshold)
        return best_score, self.get_individuals_to_genotypes_from_histogram(ordered_individuals, ordered_genotypes, genotypes_to_counts)
    
    def count_branch_and_bound_recursive(self, geometrically_sorted_individuals_and_log_likelihoods, T_nonzero, ordered_genotypes, genotypes_to_counts, number_assigned, log_mult_to_get_here, threshold):
        n = len(geometrically_sorted_individuals_and_log_likelihoods)
        number_remaining = n - number_assigned
        
        next_genotype_index = len(genotypes_to_counts)
        if next_genotype_index == len(T_nonzero):
            return log_mult_to_get_here, deepcopy(genotypes_to_counts)

        next_genotype = ordered_genotypes[next_genotype_index]
        recursive_scores_and_genotype_counts = []

        for next_number in xrange(0, number_remaining + 1):
            log_mult_for_assignment = self.log_multiplier(geometrically_sorted_individuals_and_log_likelihoods, number_assigned, T_nonzero, next_genotype, next_number, number_remaining)
            genotypes_to_counts[next_genotype] = next_number
            log_upper_bound_to_finish_from_here = self.log_upper_bound_after_assignment(geometrically_sorted_individuals_and_log_likelihoods, genotypes_to_counts, T_nonzero, number_assigned + next_number)

            new_log_mult_to_get_here = log_mult_to_get_here + log_mult_for_assignment
            if new_log_mult_to_get_here + log_upper_bound_to_finish_from_here < threshold:
                continue

            # cannot bound, so recurse
            best_remaining_score_and_genotype_counts = self.count_branch_and_bound_recursive(geometrically_sorted_individuals_and_log_likelihoods, T_nonzero, ordered_genotypes, genotypes_to_counts, number_assigned+next_number, new_log_mult_to_get_here, threshold)
            recursive_scores_and_genotype_counts.append(best_remaining_score_and_genotype_counts)

        # all outcomes of the code require returning, so the genotype
        # added at this level must be removed
        del genotypes_to_counts[next_genotype]

        if len(recursive_scores_and_genotype_counts) == 0:
            # all for this were bounded
            return -float('inf'), {}
        return max(recursive_scores_and_genotype_counts)

    # mathematical tools used by the branch and bound
    def log_multiplier(self, geometrically_sorted_individuals_and_log_likelihoods, number_assigned, T_nonzero, next_genotype, next_number, number_remaining):
        # computes the probability of drawing these individuals from
        # the distribution and producing this data
        log_prob_drawing_from_T = log_binomial_cached(number_remaining, next_number) + log_pow(T_nonzero[next_genotype], next_number)
        relevant_individuals_and_log_likelihoods = geometrically_sorted_individuals_and_log_likelihoods[number_assigned : number_assigned + next_number]
        log_like = sum( [ log_likelihoods[next_genotype] for individual, log_likelihoods in relevant_individuals_and_log_likelihoods ] )
        return log_prob_drawing_from_T + log_like

    def log_upper_bound_after_assignment(self, geometrically_sorted_individuals_and_log_likelihoods, genotypes_to_counts, T_nonzero, number_assigned):
        # compute a tight upper bound on the best remaining path
        log_best_remaining_likelihood = sum([ max( [ log_likelihoods[g] for g in T_nonzero if g not in genotypes_to_counts ] + [ -float('inf') ] ) for individual, log_likelihoods in geometrically_sorted_individuals_and_log_likelihoods[number_assigned:] ])
        # compute the sum of all probabilities not used by the
        # multinomial (take the max with 0.0 in case of very small
        # numerical error)
        remaining_prob = max(1-sum([ T_nonzero[g] for g in genotypes_to_counts]), 0.0)
        log_best_remaining_multinomial_term = log_pow( remaining_prob, len(geometrically_sorted_individuals_and_log_likelihoods) - number_assigned )
        return log_best_remaining_likelihood + log_best_remaining_multinomial_term

    def get_individuals_to_genotypes_from_histogram(self, ordered_individuals, ordered_genotypes, genotypes_to_counts):
        individuals_to_genotypes = {}
        cumulative_assigned = 0
        for genotype in ordered_genotypes:
            if genotype not in genotypes_to_counts:
                continue
            count = genotypes_to_counts[genotype]
            individuals_to_genotypes.update( dict([ (ind, genotype) for pair, ind in ordered_individuals[cumulative_assigned : cumulative_assigned + count] ]) )
            cumulative_assigned += count
        return individuals_to_genotypes

    def parameter_range_generator():
        pass
    def T_generator():
        pass
    def log_parameter_prior():
        pass

class F1Inference(PopulationInference):
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        PopulationInference.__init__(self, filename, ploidy_range, sigma_range, vars_for_MAP = ['sigma'], **Kwargs)
        self.number_of_marginalized_configs = sum( [ len([self.parameter_range_generator(ploidy = p)]) for p in self.ploidy_range ] ) - 1
    def parameter_range_generator(self, ploidy):
        for p1 in [ (i, ploidy-i) for i in xrange(0, ploidy+1) ]:
            # use symmetry to only check each pair twice by
            # constraining p2 >= p1
            for p2 in [ (i, ploidy-i) for i in xrange(p1[0], ploidy+1) ]:
                yield((p1, p2))
    def log_parameter_prior(self, ploidy, **Kwargs):
        return log(1/float(math.exp(log_binomial_cached(ploidy+2, 2))))
    def T_generator(self, ploidy, parameter):
        g_p1, g_p2 = parameter
        if ploidy % 2 == 1:
            print 'Warning: using hypergeometric distribution of offspring with an odd ploidy'
        # fill the table with every genotype of interest
        result = dict( [ ((i, ploidy-i), 0.0) for i in range(0, ploidy+1) ] )
        # ploidy should be even, so division by 2 is OK
        for x1 in xrange(0, ploidy/2+1):
            log_prob_x1 = log_binomial_cached(g_p1[0], x1) + log_binomial_cached(ploidy - g_p1[0], ploidy/2 - x1) - log_binomial_cached(ploidy, ploidy/2)
            for x2 in xrange(0, ploidy/2+1):
                log_prob_x2 = log_binomial_cached(g_p2[0], x2) + log_binomial_cached(ploidy - g_p2[0], ploidy/2 - x2) - log_binomial_cached(ploidy, ploidy/2)
                prob = math.exp(log_prob_x1 + log_prob_x2)
                offspring_geno = ( x1+x2, ploidy-(x1+x2) )
                result[ offspring_geno ] += prob
        return result

class SlowF1Inference(F1Inference):
    # overrides the branch and bound routine from ParentalInference to
    # go sloooooowwww.
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        F1Inference.__init__(self, filename, ploidy_range, sigma_range, **Kwargs)

    def count_branch_and_bound(self, individuals_to_log_likelihoods, T_nonzero, threshold, ordered_individuals, ordered_genotypes):
        return self.genotype_branch_and_bound(individuals_to_log_likelihoods, T_nonzero, threshold)

    def genotype_branch_and_bound(self, individuals_to_log_likelihoods, T_nonzero, threshold):
        ind0 = individuals_to_log_likelihoods.keys()[0]

        inds = individuals_to_log_likelihoods.keys()
        log_likes = [ individuals_to_log_likelihoods[i] for i in inds ]

        genos_to_inds = dict([ (g,i) for i,g in enumerate(T_nonzero) ])
        return self.genotype_branch_and_bound_recursive(inds, log_likes, T_nonzero, threshold, genos_to_inds, { tuple([0]*len(T_nonzero)) : (0.0, {}) }, 0)
    
    def genotype_branch_and_bound_recursive(self, index_inds, index_log_likelihoods, T_nonzero, threshold, genos_to_inds, hists_to_log_scores_and_configs, depth):
        if len(hists_to_log_scores_and_configs) == 0:
            return (-float('inf'), {})

        if depth == len(index_inds):
            return max([(log_score + self.log_hist_prob(hist, T_nonzero), config) for hist, (log_score, config) in hists_to_log_scores_and_configs.items()] )

        ind, log_likes = index_inds[depth], index_log_likelihoods[depth]
        max_remaining_after = sum( [ max(log_l.values()) for log_l in index_log_likelihoods[depth+1:] ] )

        # go through hists at this level
        new_hists_to_log_scores_and_configs = {}
        for hist, (log_score, config) in hists_to_log_scores_and_configs.items():
            for g in T_nonzero:
                i = genos_to_inds[g]
                # get new hist
                new_hist = list(hist)
                new_hist[i] += 1
                new_hist = tuple(new_hist)

                # compute the new score of this config
                new_log_score = log_score + log_likes[g]

                if new_log_score + max_remaining_after >= threshold:
                    new_config = deepcopy(config)
                    new_config[ind] = g
                    if new_hist not in new_hists_to_log_scores_and_configs:
                        new_hists_to_log_scores_and_configs[ new_hist ] = (new_log_score, new_config)
                    elif new_hist in new_hists_to_log_scores_and_configs and new_log_score > new_hists_to_log_scores_and_configs[ new_hist ][0]:
                        new_hists_to_log_scores_and_configs[ new_hist ] = (new_log_score, new_config)

        # recurse using this layer
        return self.genotype_branch_and_bound_recursive(index_inds, index_log_likelihoods, T_nonzero, threshold, genos_to_inds, new_hists_to_log_scores_and_configs, depth+1)

    def log_hist_prob(self, hist, T_nonzero):
        genos_to_counts = dict([ (g, hist[i]) for i,g in enumerate(T_nonzero) ])
        return log_probability_of_histogram_given_distribution(genos_to_counts, T_nonzero)

class all_cached_functor:
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if args in self.cache:
            return self.cache[args]
        else:
            ans = self.func(*args)
            self.cache[args] = ans
            return ans

class single_cached_functor:
    def __init__(self, func):
        self.func = func
        self.cache = {}
    def __call__(self, *args):
        if args in self.cache:
            return self.cache[args]
        else:
            if len(self.cache) > 1:
                self.cache = {}
            ans = self.func(*args)
            self.cache[args] = ans
            return ans

class F1AndParentInference(F1Inference):
    def __init__(self, parent_filename, progeny_filename, ploidy_range, sigma_range, **Kwargs):
        F1Inference.__init__(self, progeny_filename, ploidy_range, sigma_range, **Kwargs)
        self.parent_inference = SharedPloidyInference(parent_filename, ploidy_range, sigma_range, **Kwargs)
        self.number_of_marginalized_configs = sum( [ len([self.parameter_range_generator(ploidy = p)]) for p in self.ploidy_range ] ) - 1
        # note: using this type of cache may increase memory usage significantly
        self.cached_parent_inds_to_log_likelihoods = single_cached_functor(self.parent_inference.get_individuals_to_log_likelihoods)
    def parameter_range_generator(self, ploidy):
        # since the parents are labeled (unlike in the F1Inference
        # model), all (ploidy+1)*(ploidy+1) configurations must be
        # tried
        for p1 in [ (i, ploidy-i) for i in xrange(0, ploidy+1) ]:
            for p2 in [ (i, ploidy-i) for i in xrange(0, ploidy+1) ]:
                yield((p1, p2))
    def log_parameter_prior(self, ploidy, parameter, sigma):
        # sort the parent names just in case they come out in a
        # different order sometimes
        parents_to_genos = dict(zip(sorted(self.parent_inference.individuals_to_data.keys()), list(parameter)))
        # note: this would be much more efficient to use a cached functor for the parent log likelihoods
        #log_likelihood_parents = self.parent_inference.log_likelihood_of_configuration(parents_to_genos, self.parent_inference.get_individuals_to_log_likelihoods(ploidy, sigma))
        like_table = self.cached_parent_inds_to_log_likelihoods(ploidy, sigma)
        log_likelihood_parents = self.parent_inference.log_likelihood_of_configuration(parents_to_genos, like_table)
        # note: calling F1Inference log_parameter_prior is not
        # appropriate because it should return 1/(n+1)(n+1)
        #return F1Inference.log_parameter_prior(self, ploidy) + log_likelihood_parents
        return 2*log(ploidy+1) + log_likelihood_parents

class HWInference(PopulationInference):
    def __init__(self, filename, ploidy_range, sigma_range, **Kwargs):
        PopulationInference.__init__(self, filename, ploidy_range, sigma_range, vars_for_MAP = ['sigma', 'parameter'], **Kwargs)
        self.number_of_marginalized_configs = len(ploidy_range) - 1
    def parameter_range_generator(self, ploidy):
        # note: this doesn't use ploidy. it's accepted as an argument
        # because it needs to match the
        # unique_parents_for_ploidy_generator function.

        return numpy.arange(0.01, 0.99, 0.05)
    def T_generator(self, ploidy, parameter):
        x_freq = parameter
        # fill the table with every genotype of interest
        result = dict( [ ((i, ploidy-i), 0.0) for i in range(0, ploidy+1) ] )
        for g in result:
            x, y = g
            log_prob = log_binomial_cached(ploidy, x) + log_pow(x_freq, x) + log_pow(1-x_freq, ploidy-x)
            result[g] = math.exp(log_prob)
        return result
    def log_parameter_prior(self, ploidy, **Kwargs):
        return log(1 / float(ploidy+1))

