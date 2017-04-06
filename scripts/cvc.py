#!/usr/bin/env python3

import sys
import numpy as np
from flat import Flat

result = np.zeros((2, 21))
result[0,:] = np.linspace(-2,2,21)
for i, Exc in enumerate(result[0,:]):
    f = Flat(m="graphene",
             Exc=Exc,
             T=300,
             n=100,
             dt=1e-14,
             alltime=1e-8)
    d = f.run(True)
    result[1, i] = d["v_y"][0]
np.savetxt("cvc.dat", result.transpose())
