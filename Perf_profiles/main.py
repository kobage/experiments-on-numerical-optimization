#1.  To print the histogram, uncomment next 25 lines
"""
import numpy as np
import matplotlib.pyplot as plt

f, ax2 = plt.subplots(1,1)
f.set_figheight(10/2.54)

barWidth = 0.25
L_BFGS =    [ 10.676595, 11.68022,  23.676029, 14.909464,  12.161625, 21.218217,  19.737626,  16.194141, 12.861031, 15.552349 ]
cg_descent = [17.311016, 20.109558, 21.881854, 21.067816, 15.818545, 27.298047, 26.788445, 31.171577, 16.609970, 19.187344]

# Set position of bar on X axis
br1 = np.arange(len(L_BFGS))
br2 = [x + barWidth for x in br1]

# Make the plot
plt.bar(br1, L_BFGS, color='w', hatch = '', width=barWidth, edgecolor='grey', label='L-BFGS')
plt.bar(br2, cg_descent, color='w', hatch = '..',width=barWidth,  edgecolor='grey', label='cg_descent')

# Adding Xticks
ax2.set_ylabel('Fig. 3. CPU time', fontsize=10)
ax2.set_xticks([r + barWidth for r in range(len(L_BFGS))],
           ['neuron 1', 'neuron 2', 'neuron 3', 'neuron 4', 'neuron 5', 'neuron 6', 'neuron 7', 'neuron 8', 'neuron 9', 'neuron 10' ],rotation=30)
ax2.legend(prop = { "size": 7 })

plt.tight_layout()
plt.show()
"""

#2.  To print performance profiles - two profiles,  each with data from 3 files, uncomment next
"""
import numpy as np
import matplotlib.pyplot as plt

from performance_profile import build_profile

from performance_profile import build_profile_ax
a= [ 'cg_descent_6.8.txt', 'MHB with BETA = 1.txt', 'MNAG with BETA = 1.txt']
a1= ['cg_descent_6.8.txt','MHB with BETA = 1.001.txt', 'MNAG with BETA = 1.txt']

f, (ax1, ax2) = plt.subplots(2,1)
f.set_figheight(14/2.54)
f.set_figwidth(20/2.54)
plt.gcf().set_size_inches(18/2.54, 16/2.54)

#parameter #3 - Fig. number
build_profile_ax(ax1,a,10,1)
build_profile_ax(ax2,a1,50,2)

plt.tight_layout()
plt.show()
"""

#3.  To print the performance profile with data from 2 files, uncomment next
"""
from performance_profile import build_profile

from performance_profile import build_profile
a= [ 'cg_descent_6.8.txt', 'MHB with BETA = 1.txt']

build_profile(a,10,4)
"""