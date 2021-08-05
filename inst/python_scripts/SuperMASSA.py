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

from inference import *
import tempfile

def split_range_arguments(str_range, tuple_of_allowed_sizes):
    # accepts things of the form low:high:step, low:high, and val
    str_options = str_range.split(':')
    if len(str_options) not in tuple_of_allowed_sizes:
        raise Exception('this range must specify number of arguments in the following' + str(tuple_of_allowed_sizes))

    if len(str_options) == 3:
        low, high, step_size = str_options
    elif len(str_options) == 2:
        low, high = str_options
        step_size = 2
    elif len(str_options) == 1:
        low = high = str_options[0]
        step_size = 1
    # the values should be numerical, so float will handle any type of
    # result
    low, high, step_size = float(low), float(high), float(step_size)
    if low > high:
        raise Exception('Error: range requires low < high')
    return low, high, step_size

# main and profiling functions
def real_main(argv):
    try:
        opts, args = getopt.getopt(argv, '', ['inference=', 'optimal_inference', 'file=', 'ploidy_range=', 'sigma_range=', 'draw_genotypes', 'print_genotypes', 'print_genotypes_histogram', 'draw_genotype_histogram', 'heights_or_areas=', 'naive_posterior_reporting_threshold=', 'display_progress', 'save_figures=', 'f1_parent_data=', 'slow_inference', 'save_geno_prob_dist='])
    except getopt.GetoptError as error:
        print 'Command line error:' + str(error)
        sys.exit(2)

    options_map = dict( [ (flag[2:], value) for flag, value in opts ] )

    if 'file' not in options_map:
        print 'Error: you must specify a \"file\"'
        sys.exit(2)
    fname = options_map['file']

    if 'save_geno_prob_dist' not in options_map:
        print 'Error: you must specify a \"file\" to save the probability distribution of the genotypes'
        sys.exit(2)
    prob_geno_fname = options_map['save_geno_prob_dist']
        
    if 'inference' not in options_map:
        print 'Error: must specify a type of inference'
        sys.exit(2)
        
    else:
        inference = options_map['inference']
        valid_inference_set = ('ploidy', 'f1', 'hw')
        if inference not in valid_inference_set:
            print 'Error: inference must be one of these', valid_inference_set
            sys.exit(2)

    if 'ploidy_range' in options_map:
        str_ploidy_range = options_map['ploidy_range']
        low_ploidy, high_ploidy, step_size = split_range_arguments(str_ploidy_range, (1,2,3))
        low_ploidy, high_ploidy, step_size = int(low_ploidy), int(high_ploidy), int(step_size)
        ploidy_range = range( low_ploidy, high_ploidy+1, step_size )
        if len(ploidy_range) == 0:
            print 'Error: ploidy range specified is empty'
            sys.exit(2)
    else:
        ploidy_range = range(2, 16+1, 2)

    if 'sigma_range' in options_map:
        str_sigma_range = options_map['sigma_range']
        low_sigma, high_sigma, step_size = split_range_arguments(str_sigma_range, (1,3))
        low_sigma, high_sigma, step_size = float(low_sigma), float(high_sigma), float(step_size)
        # add a small constant in case the low and high are equal
        sigma_range = P.arange( low_sigma, high_sigma+1e-5, step_size )
        if len(sigma_range) == 0:
            print 'Error: sigma range specified is empty'
            sys.exit(2)
    else:
        sigma_range = [ 0.01, 0.02, 0.04, 0.08, 0.16, 0.32]

    if 'heights_or_areas' in options_map:
        if options_map['heights_or_areas'] not in ('heights','areas'):
            print 'Error: heights_or_areas must be either "heights" or "areas"'
            sys.exit(2)
        heights_or_areas = options_map['heights_or_areas']
    else:
        heights_or_areas = 'heights'

    if 'save_figures' in options_map:
        if 'draw_genotypes' not in options_map and 'draw_genotype_histogram' not in options_map:
            print 'Error: cannot save figures if there is nothing to draw'
            sys.exit(2)
        figure_prefix = options_map['save_figures']

    if 'display_progress' in options_map:
        display_progress = True
    else:
        display_progress = False

    if 'f1_parent_data' in options_map:
        if inference != 'f1':
            print 'Error: you must run f1 inference to use f1 parent data'
            sys.exit(1)

    if inference == 'f1':
        if 'f1_parent_data' in options_map:
            parent_fname = options_map['f1_parent_data']
            infer = F1AndParentInference(parent_fname, fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        else:
            if 'slow_inference' in options_map:
                infer = SlowF1Inference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
            else:
                infer = F1Inference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        if 'optimal_inference' in options_map:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal(display_progress)
        else:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal_greedy(display_progress)

        if 'naive_posterior_reporting_threshold' in options_map:
            naive_reporting_thresh = float(options_map['naive_posterior_reporting_threshold'])
            inds_to_genos = infer.naive_high_quality_genotypes(naive_reporting_thresh, inds_to_genos, parameters, prob_geno_fname)

        ploidy, sigma, parents = parameters['ploidy'], parameters['sigma'], parameters['parameter']
        if 'draw_genotype_histogram' in options_map:
            hypergeom = infer.T_generator(ploidy, parents)
            observed_histogram = histogram(inds_to_genos)
            if 'save_figures' in options_map:
                draw_histogram(observed_histogram, hypergeom, figure_prefix)
            else:
                draw_histogram(observed_histogram, hypergeom)

    elif inference == 'hw':
        infer = HWInference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        if 'optimal_inference' in options_map:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal(display_progress)
        else:
            log_joint, posterior, inds_to_genos, parameters = infer.optimal_greedy(display_progress)

        if 'naive_posterior_reporting_threshold' in options_map:
            naive_reporting_thresh = float(options_map['naive_posterior_reporting_threshold'])
            inds_to_genos = infer.naive_high_quality_genotypes(naive_reporting_thresh, inds_to_genos, parameters, prob_geno_fname)

        ploidy, sigma, allele_freq = parameters['ploidy'], parameters['sigma'], parameters['parameter']
        if 'draw_genotype_histogram' in options_map:
            hw_dist = infer.T_generator(ploidy, allele_freq)
            observed_histogram = histogram(inds_to_genos)
            if 'save_figures' in options_map:
                draw_histogram(observed_histogram, hw_dist, figure_prefix)
            else:
                draw_histogram(observed_histogram, hw_dist)

    elif inference == 'ploidy':
        if 'draw_genotype_histogram' in options_map:
            print 'Error: genotype histogram is not relevant w/o parental inference'
            sys.exit(2)

        infer = SharedPloidyInference(fname, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
        log_joint, posterior, inds_to_genos, parameters = infer.optimal(display_progress)

        if 'naive_posterior_reporting_threshold' in options_map:
            naive_reporting_thresh = float(options_map['naive_posterior_reporting_threshold'])
            inds_to_genos = infer.naive_high_quality_genotypes(naive_reporting_thresh, inds_to_genos, parameters, prob_geno_fname)

# fixme: only for CGI
    print '<h1>Results</h1>'
    print 'log Pr(D, G = g*, parameter*) =', log_joint, ', '
    print 'Pr(G = g*, parameter* | D) =', posterior, ', '
    print 'parameter*', parameters
# fixme: only for CGI
#    print ''
    print '<br>'
    if 'print_genotypes' in options_map:
        # fixme: only for CGI
        print '<br>'
        print '<h1>Predicted genotypes</h1>'
        for i,g in sorted(inds_to_genos.items()):
# fixme: only for CGI
#            print i,' ',g
            print i,'\t',g,"<br>"
# fixme: only for CGI
    print '<br>'
    if 'print_genotypes_histogram' in options_map:
        print histogram(inds_to_genos)
    if 'draw_genotypes' in options_map:
        if options_map['inference'] == 'f1' and 'f1_parent_data' in options_map:
            parent_filename = options_map['f1_parent_data']
            par_infer = SharedPloidyInference(parent_filename, ploidy_range, sigma_range, heights_or_areas = heights_or_areas)
            par_inds_to_genos = par_infer.get_inds_to_genotypes_given_parent_parameters(**parameters)
            if 'save_figures' in options_map:
                draw_genotypes(par_infer.individuals_to_data, par_inds_to_genos, 2, figure_prefix + '_parents')
            else:
                draw_genotypes(par_infer.individuals_to_data, par_inds_to_genos, 2)
            
        if 'save_figures' in options_map:
            draw_genotypes(infer.individuals_to_data, inds_to_genos, 1, figure_prefix)
        else:
            draw_genotypes(infer.individuals_to_data, inds_to_genos, 1)

    if 'save_figures' not in options_map and ( 'draw_genotypes' in options_map or 'draw_genotype_histogram' in options_map ):
        getchar()

def profile_main(argv):
    cProfile.run('real_main(' + str(argv) + ')', sort='cumulative')

def folder_prefix(s):
    return s[s.rfind('/', 0, s.rfind('/')-1)+1:]

    class cgi_tools:
        def __init__(self, reform = False):
            self.reform = reform
        def myprint(self, *args):
            if not self.reform:
                print(args)
            for i,x in enumerate(args):
                if type(x) == str:
                    self.myprint(x.replace('\n', '<BR>')),
                elif True: #x is an atom (i.e. int, float, etc.):
                    print x,
                else:
                    self.myprint(x),
                if i != len(args)-1:
                    print ' ',
            print ''
            
def cgi_main(argv):
    try:
        if '--draw_genotypes' in set(argv) or '--draw_genotype_histogram' in set(argv):
            temp_prefix='/home/mmollin/Downloads/src'
            temp, temp_fname = tempfile.mkstemp(dir=temp_prefix)
            argv.extend(['--save_figures', temp_fname])
            local_temp_fname = folder_prefix(temp_fname)
        real_main(argv)
        print '<table>'
        if '--draw_genotypes' in set(argv):
            print '<tr><td><h1>Population data</h1><td>'
            print '<img src="/' + local_temp_fname + '_scatter.png" width=600>'
            if '--f1_parent_data' in set(argv):
                print '<tr><td><h1>Parent data</h1><td>'
                print '<img src="/' + local_temp_fname + '_parents_scatter.png" width=600>'
        if '--draw_genotype_histogram' in set(argv):
            print '<tr><td><h1>Genotype frequencies</h1><td>'
            print '<img src="/' + local_temp_fname + '_hist.png" width=600>'
        print '</table>'
    except Exception as error:
        print 'CGI error: ' + str(error)
        sys.exit(2)

if __name__ == '__main__':
#    profile_main(sys.argv[1:])
    real_main(sys.argv[1:])
