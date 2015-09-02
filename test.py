
from numpy import *
from random import *
from matplotlib.pyplot import *
rcdefaults()

p1 = plot([1,2,3])
p2 = plot([3,2,1])
p3 = plot([2,3,1])
legend([p2, p1], ["line 2", "line 1"])
show()