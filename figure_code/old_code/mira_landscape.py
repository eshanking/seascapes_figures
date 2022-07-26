# cefotaxime landscape from mira at al
# Mira PM, Crona K, Greene D, Meza JC, Sturmfels B, Barlow M. 
# Rational Design of Antibiotic Treatment Plans: A Treatment Strategy for Managing Evolution and Reversing Resistance. 
# Planet PJ, ed. PLoS ONE. 2015;10(5):e0122283. doi:10.1371/journal.pone.0122283

# cefotaxime concentration = 0.05 ug/mL

import matplotlib.pyplot as plt
import numpy as np
from seascapes_figures.classes.population_class import Population

p = Population()

cef_landscape = [.160, 0.185, 1.653, 
                1.936, 0.085, 0.225, 
                1.969, 0.140, 2.295, 
                0.138, 2.348, 0.119, 
                0.092, 0.203, 2.269, 2.412]

cef_landscape = np.array(cef_landscape)

fig,ax = plt.subplots()

ax = p.plot_landscape(fitness=cef_landscape)