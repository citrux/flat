#!/usr/bin/env python3

import numpy as np
from flat import Flat

freqs = np.linspace(5.25e14, 5.35e14, 3)
waves = lambda freq: [[1, 0, 0, 1e12, 0], [10, 0, 0, freq, 0]]
population = []

for freq in freqs:
    f = Flat(m="bigraphene_0.01_2", waves=waves(freq), T=300, n=100, dt=1e-15, alltime=1e-9)
    d = f.run()
    population.append(d["population"])

population = np.array(population)
data = np.hstack([freqs.reshape(len(freqs), 1), population[:,:,0]])
np.savetxt("ps.dat", data)
