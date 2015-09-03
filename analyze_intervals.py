
from scipy.stats import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

with open('backward_intervals_2.txt', 'rb') as f:
    intervals_info = pickle.load(f)
    
gints = [i[4] for i in intervals_info] #get all generation times
#n, bins, patches = plt.hist(gints, 50, normed=1, histtype='stepfilled')
#plt.show()

intrinsic_gmean = 3*1/.5 + 3*1/.5
print 'found mean', np.mean(gints)
print 'real mean', intrinsic_gmean
