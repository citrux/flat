#!/usr/bin/env python3

import re
import sys
import numpy as np
from subprocess import Popen, PIPE
#from itertools import product
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


iteration = 0
program = "./flat"
if sys.platform.startswith("win32"):
    program = "flat.exe"

def parse(name, text):
    mean, std = re.search(name + r"\s*=\s*(\S*)\s*\+/-\s*(\S*)", text).groups()
    return [float(mean), float(std)]

def calculate(prefix, mode, Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau):
    global iteration
    iteration += 1
    print("%d" % iteration)
    args = "%e %e %e %e %e %e %e %e %d %e %e\n" % (Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime)
    out, err = Popen([program], shell=True, stdin=PIPE, stdout=PIPE).communicate(input=args.encode("utf8"))
    data = out.decode("utf8")
    print(data)
    return {
        "v_x": parse("v_x", data),
        "v_y": parse("v_y", data),
        "power": parse("power", data),
        "tau": parse("tau", data),
        # "upper_band": parse("upper", data),
        # "lower_band": parse("lower", data)
    }

def cvc():
    Exc     = 0
    Eyc     = 0
    H       = 0
    Ex      = 0
    Ey      = 0
    omega   = 0
    phi     = 0
    T       = 300
    n       = 16
    dt      = 5e-14
    alltime = 4e-9
    tau     = 3e-12
    Excl = np.linspace(0, 1, 21)
    one = [calculate("CVC1", "one_band", Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for Exc in Excl]
    vx = np.array([x["v_x"] for x in one])
    vy = np.array([x["v_y"] for x in one])
    plt.plot(Excl, vx[:,0])
    plt.plot(Excl, vy[:,0])
    plt.show()

def freq():
    Exc     = 1
    Eyc     = 0
    H       = 0
    Ex      = 0.3
    Ey      = 0
    omegal  = np.linspace(5e9, 5e11, 100)
    phi     = np.pi/2
    T       = 300
    n       = 16
    dt      = 1e-13
    alltime = 1e-8
    tau     = 3e-12
    # однозонное приближение
    zer = [calculate("CVC1", "one_band", Exc, Eyc,   0, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
    one = [calculate("CVC1", "one_band", Exc, Eyc, 200, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
    two = [calculate("CVC1", "one_band", Exc, Eyc, 400, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
    # двухзонное
    # two_bands = [calculate("CVC2", "two_bands", Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
    #print(one_band)
    #print(two_bands)
    power0 = np.array([x["power"] for x in zer])
    power1 = np.array([x["power"] for x in one])
    power2 = np.array([x["power"] for x in two])
    # power2 = np.array([x["power"] for x in two_bands]) * 10e12 * 1.6e-19 / dt
    # tau1 = np.array([x["tau"] for x in one_band])
    # tau2 = np.array([x["tau"] for x in two_bands])
    # lower = np.array([x["lower_band"] for x in two_bands])
    # upper = np.array([x["upper_band"] for x in two_bands])
    # plt.errorbar(omegal, power1[:,0], yerr=power1[:,1], label="Exc=%.2f" % Exc)
    plt.plot(omegal, power0[:,0], label="H = 0")
    plt.plot(omegal, power1[:,0], label="H = 200")
    plt.plot(omegal, power2[:,0], label="H = 400")
    # plt.errorbar(omegal, power2[:,0], yerr=power2[:,1], label="two bands")
    plt.grid()
    plt.legend(loc='upper right')
    plt.title("Absorption")
    plt.show()

if __name__ == '__main__':
    freq()
    #cvc()
