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

import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')

import numpy
import pylab as P
import sys

def draw_histogram(experimental_geno_hist, theoretical_geno_hist, save_file_prefix = '', file_type = 'png'):
    if len(experimental_geno_hist) == 0:
        raise Exception('cannot draw histogram of no genotypes')

    if save_file_prefix == '':
        P.ion()
        P.figure(0)
    else:
        P.figure(0)
    P.clf()

    total = sum( experimental_geno_hist.values() )

    P.bar( *zip(* [ (g[0], theoretical_geno_hist[g]) for g in sorted(theoretical_geno_hist)] ), width = 1, color = 'blue', alpha = 0.9 , label = 'Expected')
    P.bar( *zip(* [ (g[0], experimental_geno_hist[g] / float(total)) for g in sorted(experimental_geno_hist)] ), width = 1, color = 'red', alpha = 0.7, label = 'Observed' )

    P.xlabel('Doses of allele 1', size='x-large')
    P.ylabel('Frequencies', size='x-large')

    P.legend(loc='upper right')
    
    if save_file_prefix == '':
        P.draw()
    else:
        P.savefig(save_file_prefix + '_hist.' + file_type, dpi=600)

def draw_genotypes(individuals_to_data, individuals_to_genotypes, figure_number, save_file_prefix = '', file_type = 'png'):
    # if set(individuals_to_data) is a superset of
    # set(individuals_to_genotypes), then only draw the points for the
    # assigned genotypes
    if len(individuals_to_genotypes) == 0:
        raise Exception('cannot draw scatterplot of no genotypes')

    base_colors = ['red', 'blue', 'green', 'black', 'pink', 'orange', 'purple', 'yellow', 'cyan', 'magenta', 'gray']
    base_markers = [ 's' , 'o' , '^', '+' , 'd' , 'h' , '>' , 'p' , '<' , 'v' , 'x', '1' , '2' , '3' , '4', 'None' , ' ' , '' , '$...$']

    # make base colors and markers the same length
    base_display_number = min(len(base_colors), len(base_markers))
    base_colors = base_colors[:base_display_number]
    base_markers = base_markers[:base_display_number]

    # pair every color with every marker
    colors = []
    markers = []
    for j in xrange(0, base_display_number):
        for i,c in enumerate(base_colors):
            colors.append(c)
            markers.append(base_markers[(j+i)%base_display_number])

    if save_file_prefix == '':
        P.ion()
        P.figure(figure_number)
    else:
        P.figure(figure_number)
    P.clf()

    if len(set(individuals_to_genotypes.values())) > min(len(colors), len(markers)):
        print 'Error: cannot plot more than', min(len(colors), len(markers)), 'distinct colors at this time (automated color generation could lower contrast between adjacent classes)'
        sys.exit(2)
    
    genos_to_x_and_ys = {}
    for individual, geno in individuals_to_genotypes.items():
        if geno not in genos_to_x_and_ys:            genos_to_x_and_ys[geno] = individuals_to_data[individual]
        else:
            genos_to_x_and_ys[geno].extend(individuals_to_data[individual])

    x_max = max([ max([x_and_y[0] for x_and_y in genos_to_x_and_ys[g]]) for g in genos_to_x_and_ys ])
    y_max = max([ max([x_and_y[1] for x_and_y in genos_to_x_and_ys[g]]) for g in genos_to_x_and_ys ])
    x_or_y_max = max(x_max, y_max)
    for g, c in zip(sorted(genos_to_x_and_ys), colors):
        # make the theoretical series for that genotype
        if g[0] > 0.0:
            r_limit_x = x_max / g[0]
        else:
            r_limit_x = 1.0e10

        if g[1] > 0.0:
            r_limit_y = y_max / g[1]
        else:
            r_limit_y = 1.0e10

        r_limit = min( r_limit_x, r_limit_y )
        s = [ (0,0), (g[0]*r_limit, g[1]*r_limit) ]
        P.plot( *zip(*s),  c = c, linewidth=2)

    for g, c, m in zip(sorted(genos_to_x_and_ys), colors, markers):
        P.scatter( *zip(*genos_to_x_and_ys[g]), c = c , marker = m, s=40)

    # this adds a little buffer for the figure; note that pylab.axes should probably be used somehow
    P.scatter([-0.05*x_max,-0.05*x_max,x_max,x_max],[-0.05*y_max,y_max,-0.05*y_max,y_max], visible = False)
    P.xlabel('Intensity of allele 1', size='x-large')
    P.ylabel('Intensity of allele 2', size='x-large')
    P.xlim(0,x_or_y_max)
    P.ylim(0,x_or_y_max)

    if save_file_prefix == '':
        P.draw()
    else:
        P.savefig(save_file_prefix + '_scatter.' + file_type, dpi=600)

def getchar():
    print ''
    Var = raw_input("press a key")

def genotype_histogram(P, genotypes):
    result = dict( [ ((i, P-i), 0.0) for i in xrange(0, P+1) ] )
    
    for g in genotypes:
        result[g] += 1
    return result

def histogram(dictionary):
    # builds a histogram of the outputs of the dictionary
    result = {}
    for val in dictionary.values():
        if val in result:
            result[val] += 1
        else:
            result[val] = 1
    return result

def nonzero_distribution(dist):
    return dict([(key, value) for key, value in dist.items() if value != 0.0])

def indent(i_val):
    for i in range(0, i_val):
        print '  ',
