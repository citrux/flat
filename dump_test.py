#!/usr/bin/env python3
from flat import Flat
f = Flat(dumping=True,
         m="bigraphene",
         Exc=1,
         Hc=200,
         Ex=0.3,
         omega=5e10,
         T=300,
         n=16,
         dt=1e-13,
         alltime=1e-8)
f.run(True)
