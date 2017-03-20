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

def magnetic():
    Exc     = 1
    Eyc     = 0
    Hl      = np.linspace(0, 1000, 41)
    Ex      = 0.1
    Ey      = 0
    omega   = 1e11
    phi     = 0
    T       = 300
    n       = 1000
    dt      = 5e-14
    alltime = 1e-9
    tau     = 3e-12
    # однозонное приближение
    one_band = [calculate("CVC1", "one_band", Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for H in Hl]
    # двухзонное
    two_bands = [calculate("CVC2", "two_bands", Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for H in Hl]
    #print(one_band)
    #print(two_bands)
    power1 = np.array([x["power"] for x in one_band]) * 10e12 * 1.6e-19 / dt
    power2 = np.array([x["power"] for x in two_bands]) * 10e12 * 1.6e-19 / dt
    tau1 = np.array([x["tau"] for x in one_band])
    tau2 = np.array([x["tau"] for x in two_bands])
    lower = np.array([x["lower_band"] for x in two_bands])
    upper = np.array([x["upper_band"] for x in two_bands])
    plt.errorbar(Hl, power1[:,0], yerr=power1[:,1], label="one band")
    plt.errorbar(Hl, power2[:,0], yerr=power2[:,1], label="two bands")
    plt.grid()
    plt.legend()
    plt.title("Absorption")
    plt.show()

    plt.errorbar(Hl, tau1[:,0], yerr=tau1[:,1], label="one band")
    plt.errorbar(Hl, tau2[:,0], yerr=tau2[:,1], label="two bands")
    plt.grid()
    plt.legend()
    plt.title("Average relaxation time")
    plt.show()

    plt.errorbar(Hl, lower[:,0], yerr=lower[:,1], label="lower")
    plt.errorbar(Hl, upper[:,0], yerr=upper[:,1], label="upper")
    plt.grid()
    plt.legend()
    plt.title("Time in bands")
    plt.show()

def freq():
    Exc     = 0
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
    for Exc in [0]: #np.linspace(0, 2, 5):
        # однозонное приближение
        one = [calculate("CVC1", "one_band", Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
        two = [calculate("CVC1", "one_band", Exc, Eyc, 500, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
        # двухзонное
        # two_bands = [calculate("CVC2", "two_bands", Exc, Eyc, H, Ex, Ey, omega, phi, T, n, dt, alltime, tau) for omega in omegal]
        #print(one_band)
        #print(two_bands)
        power1 = np.array([x["power"] for x in one])
        power2 = np.array([x["power"] for x in two])
        # power2 = np.array([x["power"] for x in two_bands]) * 10e12 * 1.6e-19 / dt
        # tau1 = np.array([x["tau"] for x in one_band])
        # tau2 = np.array([x["tau"] for x in two_bands])
        # lower = np.array([x["lower_band"] for x in two_bands])
        # upper = np.array([x["upper_band"] for x in two_bands])
        # plt.errorbar(omegal, power1[:,0], yerr=power1[:,1], label="Exc=%.2f" % Exc)
        plt.plot(omegal, power1[:,0], label="without H")
        plt.plot(omegal, power2[:,0], label="with H")
        # plt.errorbar(omegal, power2[:,0], yerr=power2[:,1], label="two bands")
    plt.grid()
    plt.legend(loc='upper right')
    plt.title("Absorption")

    plt.show()

    # plt.errorbar(omegal, tau1[:,0], yerr=tau1[:,1], label="one band")
    # plt.errorbar(omegal, tau2[:,0], yerr=tau2[:,1], label="two bands")
    # plt.grid()
    # plt.legend()
    # plt.title("Average relaxation time")
    # plt.show()

    # plt.errorbar(omegal, lower[:,0], yerr=lower[:,1], label="lower")
    # plt.errorbar(omegal, upper[:,0], yerr=upper[:,1], label="upper")
    # plt.grid()
    # plt.legend()
    # plt.title("Time in bands")
    # plt.show()


def distrib():
    Exc     = 0
    Eyc     = 0
    H       = 500
    Ex      = 0.1
    Ey      = 0
    omega1   = 1.5e11
    omega2   = 5e11
    phi     = 0
    T       = 300
    n       = 100
    dt      = 1e-13
    alltime = 1e-8
    tau     = 3e-12
    dtype = np.dtype([("px", np.float32), ("py", np.float32)])
    c1 = calculate("resonance", "one_band", Exc, Eyc, H, Ex, Ey, omega1, phi, T, n, dt, alltime, tau)
    data1 = np.fromfile("resonance.bin", dtype=dtype)
    c2 = calculate("noresonance", "one_band", Exc, Eyc, H, Ex, Ey, omega2, phi, T, n, dt, alltime, tau)
    data2 = np.fromfile("noresonance.bin", dtype=dtype)

    plt.subplot(1,2,1)
    plt.hist2d(data1["px"], data1["py"], bins=400, norm=LogNorm(), cmap=plt.get_cmap("viridis"))
    plt.title("resonance")
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.subplot(1,2,2)
    plt.hist2d(data2["px"], data2["py"], bins=400, norm=LogNorm(), cmap=plt.get_cmap("viridis"))
    plt.title("no resonance")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()

def distrib2():
    Exc     = 0.5
    Eyc     = 0
    H       = 500
    Ex      = 0.1
    Ey      = 0.1
    omega1   = 1.2e11
    omega2   = 2.5e11
    omega3   = 5e11
    phi     = np.pi/2
    T       = 300
    n       = 100
    dt      = 1e-13
    alltime = 4e-8
    tau     = 3e-121
    dtype = np.dtype([("px", np.float32), ("py", np.float32)])
    c1 = calculate("resonance", "one_band", Exc, Eyc, H, Ex, Ey, omega1, phi, T, n, dt, alltime, tau)
    c2 = calculate("resonance2", "one_band", Exc, Eyc, H, Ex, Ey, omega2, phi, T, n, dt, alltime, tau)
    c3 = calculate("noresonance", "one_band", Exc, Eyc, H, Ex, Ey, omega2, phi, T, n, dt, alltime, tau)

    plt.subplot(1,3,1)
    data = np.fromfile("resonance.bin", dtype=dtype)
    plt.hist2d(data["px"], data["py"], bins=400, norm=LogNorm(), cmap=plt.get_cmap("viridis"))
    plt.title("resonance")
    xlim = plt.xlim()
    ylim = plt.ylim()
    plt.subplot(1,3,2)
    data = np.fromfile("resonance2.bin", dtype=dtype)
    plt.hist2d(data["px"], data["py"], bins=400, norm=LogNorm(), cmap=plt.get_cmap("viridis"))
    plt.title("resonance2")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.subplot(1,3,3)
    data = np.fromfile("noresonance.bin", dtype=dtype)
    plt.hist2d(data["px"], data["py"], bins=400, norm=LogNorm(), cmap=plt.get_cmap("viridis"))
    plt.title("no resonance")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()

if __name__ == '__main__':
    freq()
    # distrib2()
