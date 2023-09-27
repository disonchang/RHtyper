#!/usr/bin/env python

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys, logging
from scipy import stats
#from sklearn.metrics import mean_squared_error

log = logging.getLogger()
logging.basicConfig(level=logging.WARN)
np.random.seed(1)


def max_tstat(x):
    '''Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1'''
    sys.setrecursionlimit(10000)
    ### determine the difference from the median 
    x0 = x - np.mean(x)

    m = len(x0)
    y = np.cumsum(x0)
    ### use cumulative difference to determine the maximum boundary of change points (index of the position)
    b0, b1 = np.argmin(y), np.argmax(y)
    ### switch the index and make the index in order from low to high
    i, j = min(b0, b1), max(b0, b1)

    ### if two values equal ... added 0117_2022
    #if b0==b1: j=i

    ### 
    Y=x[i+1:j]
    Z=np.delete(x, np.arange(i+1,j))

    max_t, pval=stats.ttest_ind(Y,Z, equal_var=True, nan_policy='omit')
    max_start=i+1
    max_end=j
    return abs(max_t), max_start, max_end


def cbs(x, shuffles=1000, p=0.05, w=10):
    '''Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.'''
    max_t, max_start, max_end=max_tstat(x)
    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end, 1
    if np.isnan(max_t):
        return False, max_t, max_start, max_end, 1
    #if max_end==max_start:   ###... added 0117_2022
    #    return False, max_t, max_start, max_end, 1

    ### adjust terminal block
    if max_start < w:
        max_start = 0
    if len(x)-max_end < w:
        max_end = len(x)

    lowstat_count=0
    alpha = shuffles*p
    xt = x.copy()
    
    for i in range(shuffles):
        np.random.shuffle(xt)
        s_max_t, s_max_startIx, s_max_endIx = max_tstat(xt)
        
        #if s_max_t is None: ####... added 0117_2022
        #   continue

        if s_max_t > max_t:
           lowstat_count += 1
           
        if lowstat_count > alpha:
            ### number of inconsistent call is larger than FDR cutoff
            perm_p=1-lowstat_count/float(shuffles)
            return False, max_t, max_start, max_end, perm_p
    perm_p=1-lowstat_count/float(shuffles)
    return True, max_t, max_start, max_end, perm_p



def tstat(x, i):
    '''Return the segmentation statistic t testing if i is a (one-sided)  breakpoint in x'''
    n = len(x)
    s0 = x[:i]
    s1 = x[i:]
    t,p =stats.ttest_ind(s0,s1, equal_var=True, nan_policy='omit')
    return abs(t)

def tstat2(x,i,j):
    n = len(x)
    s0 = x[i+1:j]
    s1 = np.delete(x, np.arange(i+1,j))
    t,p =stats.ttest_ind(s0,s1, equal_var=True, nan_policy='omit')
    return abs(t)

def recursive_segment(x, start, end, L=[], shuffles=1000, p=.05, w=10):
    '''
        Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    '''
    sigseg_found, max_t, max_s, max_e, perm_p = cbs(x[start:end], shuffles=shuffles, p=p, w=w)
    #print '[rsegment]', 'Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(start, end, start+max_s, start+max_e, max_t, sigseg_found)
    #log.info('Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(start, end, start+max_s, start+max_e, max_t, sigseg_found))
    ### 3 chosen here as one window is 100bp with step size 50, 3 window will be 150bp
    if not sigseg_found or (max_e-max_s) < w or (max_e-max_s) == (end-start):
    #if not sigseg_found or (max_e-max_s) == (end-start):
         #print 'append', start, end
         L.append((start, end))
    else:
        if max_s > 0:
            ### do segmentation on the left block
            #print 'left', start, start+max_s-1
            recursive_segment(x, start, start+max_s-1, L, w=w)
        if max_e-max_s > 0:
            ### do segmentation on the middle block
            #print 'middle', start+max_s, start+max_e, 
            recursive_segment(x, start+max_s, start+max_e, L, w=w)
        if start+max_e < end:
            #print 'center', start+max_e, end
            ### do segmentation on the central block
            recursive_segment(x, start+max_e, end, L, w=w)
    return L




def segment(x, shuffles=10, p=.05, w=10):
    '''Segment the array x, using significance test based on shuffles rearrangements and significance level p
    '''
    start = 0 ### start index 
    end = len(x)-1 ### end index
    L = []
    recursive_segment(x, start, end, L, shuffles=shuffles, p=p, w=w)
    return L


def validate(x, L, shuffles=1000, p=.01):
    S = [r[0] for r in L]+[len(x)]
    SV = [0]
    left = 0
    for test, s in enumerate(S[1:-1]):
        t = tstat(x[S[left]:S[test+2]], S[test+1]-S[left])
        if t is None: continue
        #log.info('Testing validity of {} in interval from {} to {} yields statistic {}'.format(S[test+1], S[left], S[test+2], t))
        threshold = 0
        thresh_count = 0
        site = S[test+1]-S[left]
        xt = x[S[left]:S[test+2]].copy()
        flag = True
        for k in range(shuffles):
            np.random.shuffle(xt)
            threshold = tstat(xt, site)
            if threshold > t:
                thresh_count += 1
            if thresh_count >= p*shuffles:
                flag = False
                log.info('Breakpoint {} rejected'.format(S[test+1]))
                break
        if flag:
            #log.info('Breakpoint {} accepted'.format(S[test+1]))
            SV.append(S[test+1])
            left += 1
    SV.append(S[-1])
    return SV


def generate_normal_time_series(num, minl=50, maxl=1000):
    '''Generate a time series with num segments of minimal length minl and maximal length maxl.  Within a segment,
    data is normal with randomly chosen, normally distributed mean between -10 and 10, variance between 0 and 1.
    '''
    data = np.array([], dtype=np.float64)
    partition = np.random.randint(minl, maxl, num)
    for p in partition:
        mean = np.random.randn()*10
        var = np.random.randn()*1
        if var < 0:
            var = var * -1
        tdata = np.random.normal(mean, var, p)
        data = np.concatenate((data, tdata))
    return data


def draw_segmented_data(data, S, title=None):
    '''Draw a scatterplot of the data with vertical lines at segment boundaries and horizontal lines at means of 
    the segments. S is a list of segment boundaries.'''
    j=sns.scatterplot(range(len(data)),data,color='black',size=.1,legend=None)
    for x in S:
        j.axvline(x)
    for i in range(1,len(S)):
        j.hlines(np.mean(data[S[i-1]:S[i]]),S[i-1],S[i],color='green')
    j.set_title(title)
    j.get_figure().set_size_inches(16,4)
    return j

if __name__ == '__main__':

    log.setLevel(logging.INFO)
    sample = generate_normal_time_series(3, minl=5, maxl=10)
    L = segment(sample, shuffles=10)
    S = validate(sample, L)
    ax = draw_segmented_data(sample,  S, title='Circular Binary Segmentation of Data')
    ax.get_figure().savefig('plot.png')
