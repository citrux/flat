#!/usr/bin/env python3

import sys
import numpy as np
from matplotlib import pyplot as plt
from flat import Flat

Excl    = np.linspace(0, 2, 6)
Hcl     = np.linspace(20, 5000, 50)
for Exc in Excl:
    line = []
    for Hc in Hcl:
        f = Flat(m="bigraphene_0.01_1",
                 Exc=Exc,
                 Hc=Hc,
                 waves=[[1, 0, 0, 1e12, 0]],
                 T=300,
                 n=100,
                 dt=1e-13,
                 alltime=1e-9)
        line.append(f.run(True)["power"][1][0])
    plt.plot(Hcl, line, label="Exc=%.2f" % Exc)
plt.legend()
plt.gca().set_ylim(bottom=0)
plt.savefig("magnetic.png")

