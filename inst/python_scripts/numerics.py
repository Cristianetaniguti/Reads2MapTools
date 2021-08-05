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

import math
import random
import sys

def log(val_f):
    if val_f == 0:
        return -float('inf')
    return math.log(val_f)

def log_pow(base, exponent):
    if exponent == 0:
        return 0
    return exponent*log(base)

def normal_pdf(deviation, sigma):
    return 1/math.sqrt(2*math.pi*sigma)*math.exp(- deviation*deviation / (2*sigma*sigma))

def log_normal_pdf(deviation, sigma):
    return -log(math.sqrt(2*math.pi*sigma)) + -deviation*deviation / (2*sigma*sigma)

#def log_squared_normal_pdf(deviation, sigma):
#    return -2*log(math.sqrt(2*math.pi*sigma)) + -2*deviation*deviation / (2*sigma*sigma)

# fixme: some error here
# def binomial(n,x):
#     if n > x:
#         return 0
#     return prod( [ float(i) for i in xrange(x+1, n+1) ] ) / prod([float(i) for i in xrange(2, n-x+1)])

def log_binomial(n, x):
    if x > n:
        return log(0)
    return sum( [ log(i) for i in xrange(x+1, n+1) ] ) - sum([log(i) for i in xrange(2, n-x+1)])

log_sums_cache = {}
tot = 0.0
for i in xrange(1, 400):
    tot += log(i)
    log_sums_cache[i] = tot
# not accurate for log sums, but for factorial interpretation
log_sums_cache[0] = 0 
#print log_sums_cache

def log_binomial_cached(n, x):
    if n in log_sums_cache:
        if x > n:
            return log(0)
        else:
            # if in the cache, return the value
#            print 'all vals in cache', n, x
#            print log_sums_cache[n] - log_sums_cache[x] - log_sums_cache[n-x]
            return log_sums_cache[n] - log_sums_cache[x] - log_sums_cache[n-x]
    else:
        print 'fixme: implement cache extension'
        print 'wanted log_binomial of', n, x
        raise Exception('fixme: remove this later')
        sys.exit(1)
        # extend the cache
        for i in xrange(max(), max(n,x)+1):
            tot += log(i)
            log_sums_cache[i] = tot

def log_add(logA,logB):
    if logA == log(0):
        return logB
    if logA<logB:
        return log_add(logB,logA)
    return log( 1 + math.exp(logB-logA) ) + logA

def log_sum(lst):
    # given a list of floats [log(a), log(b), log(c), ...]  compute
    # log(a+b+c+d+...)
    
    # note: this could be faster by not calling log_add

    result = -float('inf')
    for i in lst:
        result = log_add(i, result)
    return result

def sample_from_distribution(dist):
    u = random.uniform(0.0, 1.0)
    
    cum = 0.0
    for outcome in dist:
        cum += dist[outcome]
        if cum >= u:
            return outcome
    print 'warning: sum of cdf from sample_from_distribution was < 1.0'

def sample_from_normal(mean, var):
    return random.normalvariate(mean, var)

def log_probability_of_histogram_given_distribution(counts_histogram, probability_dist):
    # this uses a log multinomial
    total_counts = sum( counts_histogram.values() )
    
    result = 0.0
    total_remaining = total_counts
    for g,c in counts_histogram.items():
        if c == 0:
            continue
        result += log_binomial_cached( total_remaining, c )
        if g not in probability_dist:
            return -float('inf')
        result += log_pow( probability_dist[g], c)
        total_remaining -= c
    return result

