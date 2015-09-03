
from scipy.stats import *
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(gamma.ppf(0.01, a),gamma.ppf(0.99, a), 100)
print gamma.ppf(0.01, a)
plt.plot(x, gamma.pdf(x, a),'r-', lw=5, alpha=0.6, label='gamma pdf')

import pickle

with open('backward_intervals.txt', 'rb') as f:
    intervals_info = pickle.load(f)
    
gints = [i[4] for i in intervals_info]
n, bins, patches = plt.hist(gints, 50, normed=1, histtype='stepfilled')
plt.show()