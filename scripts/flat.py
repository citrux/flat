import sys
import re
import json
from subprocess import Popen, PIPE

class Flat:
    def __init__(self, dumping=False, m="",
                       Exc=0, Eyc=0, Hc=0,
                       waves=[], T=0,
                       n=0, dt=0, alltime=0):
        self.program = "./flat"
        if sys.platform.startswith("win32"):
            self.program = "./flat.exe"
        if dumping:
            self.program += " -d"
        self.args = m
        self.args += "\n%e %e %e" % (Exc, Eyc, Hc)
        self.args += "\n%d" % (len(waves))
        for Ex, Ey, H, omega, phi in waves:
            self.args += "\n%e %e %e %e %e" % (Ex, Ey, H, omega, phi)
        self.args += "\n%e %d %e %e\n" % (T, n, dt, alltime)

    def run(self, verbose=False):
        out, err = Popen([self.program], shell=True, stdin=PIPE, stdout=PIPE)\
                   .communicate(input=self.args.encode("utf8"))
        data = out.decode("utf8")
        if verbose:
            print(data)
        return json.loads(data)
