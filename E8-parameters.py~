from scipy.stats import chi2
from math import sqrt, log, ceil

def ErrorRate(n, l, q, g, sigma, t):
    div_t = 0
    avg_t = 0.
    for i in xrange(-2**(t-1), 2**(t-1) + 1):
        div_t += (1.*i - 0)**2
    div_t /= 2.**t + 1
    #print t, div_t
    s = sqrt(n*l*sigma**2 * (sigma**2 + div_t) + n*l*sigma**4 + sigma**2)
    dis = (q - 1.) / 2 - sqrt(2.)*(q*1. / g + 1.)
    #print dis, s
    pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(n/16., 2)
    return pr

def bandwidth(n, l, q, g, sigma, t):
    seedLen = 256.
    #print (seedLen + l*n*ceil(log(q/2.**t, 2))) / 8., (l*n*ceil(log(q/2.**t, 2)) + ceil(log(g, 2))*n) / 8.
    return (seedLen + 1.*n*l*ceil(log(q/2.**t, 2)) + 1.*n*l*ceil(log(1.*q, 2)) + ceil((log(g, 2)))*n)/8.

def keySize(n, l, q, g, sigma, t):
    return n*1./2.

n, l, q, g, sigma, t = 1024, 1, 12289, 2**3, sqrt(4.), 3
print "Err rate = ", ErrorRate(n, l, q, g, sigma, t)
print "bandwidth = ", bandwidth(n, l, q, g, sigma, t)
print "keySize = ", keySize(n, l, q, g, sigma, t)
