#!/usr/bin/env python3
import numpy as np
from matplotlib import pyplot as plt
from flat import Flat
f = Flat(dumping=True,
         m="graphene",
         Exc=1,
         T=300,
         n=100,
         dt=1e-13,
         alltime=1e-7)
f.run(True)
thetas = np.loadtxt("theta.dat")
plt.hist(thetas, 40)
plt.show()
