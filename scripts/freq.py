#!/usr/bin/env python3

import sys
import numpy as np
from flat import Flat

n_e = 3
n_o = 10
Excl = [0, 0.2, 0.5]
omegal = np.linspace(5e11/n_o, 5e11, n_o)
woH = np.zeros((n_e, n_o))
woHt = np.zeros((n_e, n_o))
for i, Exc in enumerate(Excl):
    for j, omega in enumerate(omegal):
        f = Flat(m="graphene",
                 Exc=Exc,
                 Hc=0,
                 Ex=0.3,
                 omega=omega,
                 T=100,
                 n=100,
                 dt=1e-14,
                 alltime=1e-8)
        d = f.run(True)
        woH[i, j] = d["power"][0]
        woHt[i, j] = d["tau"][0]
wH = np.zeros((n_e, n_o))
'''
for i, Exc in enumerate(Excl):
    for j, omega in enumerate(omegal):
        f = Flat(m="graphene",
                 Exc=Exc,
                 Hc=400,
                 Ex=0.3,
                 omega=omega,
                 T=300,
                 n=400,
                 dt=1e-13,
                 alltime=1e-8)
        wH[i, j] = f.run(True)["power"][0]
'''
np.savetxt("woH.dat", np.hstack([omegal.reshape(n_o, 1), woH.transpose()]), header="# omega\t" + "\t".join("Exc=%.2f" % Exc for Exc in Excl))
np.savetxt("woHt.dat", np.hstack([omegal.reshape(n_o, 1), woHt.transpose()]), header="# omega\t" + "\t".join("Exc=%.2f" % Exc for Exc in Excl))
#np.savetxt("wH.dat", np.hstack([omegal.reshape(n_o, 1), wH.transpose()]), header="# omega\t" + "\t".join("Exc=%.2f" % Exc for Exc in Excl))

'''
for i, Exc in enumerate(Excl):
    for j, omega in enumerate(omegal):
        f = Flat(m="bigraphene",
                 Exc=Exc,
                 Hc=0,
                 Ex=0.1,
                 omega=omega,
                 T=300,
                 n=400,
                 dt=1e-13,
                 alltime=1e-7)
        woH[i, j] = f.run(True)["power"][0]
for i, Exc in enumerate(Excl):
    for j, omega in enumerate(omegal):
        f = Flat(m="bigraphene",
                 Exc=Exc,
                 Hc=500,
                 Ex=0.1,
                 omega=omega,
                 T=300,
                 n=400,
                 dt=1e-13,
                 alltime=1e-8)
        wH[i, j] = f.run(True)["power"][0]
'''
#np.savetxt("biwoH.dat", np.hstack([omegal.reshape(n_o, 1), woH.transpose()]), header="# omega\t" + "\t".join("Exc=%.2f" % Exc for Exc in Excl))
#np.savetxt("biwH.dat", np.hstack([omegal.reshape(n_o, 1), wH.transpose()]), header="# omega\t" + "\t".join("Exc=%.2f" % Exc for Exc in Excl))
