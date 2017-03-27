#!/usr/bin/env python3

import sys
import numpy as np
from matplotlib import pyplot as plt
from flat import Flat

Excl    = np.linspace(0, 2, 6)
Hcl     = np.linspace(20, 1000, 50)
for Exc in Excl:
    line = []
    for Hc in Hcl:
        f = Flat(m="bigraphene",
                 Exc=Exc,
                 Hc=Hc,
                 Ex=0.3,
                 omega=5e10,
                 T=300,
                 n=1000,
                 dt=1e-13,
                 alltime=1e-8)
        line.append(f.run(True)["power"][0])
    plt.plot(Hcl, line)
plt.legend(loc="upper right")
plt.gca().set_ylim(bottom=0)
plt.savefig("magnetic.png")

