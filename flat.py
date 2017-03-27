import sys
import re
from subprocess import Popen, PIPE

class Flat:
    def __init__(self, dumping=False, m="",
                       Exc=0, Eyc=0, Hc=0,
                       Ex=0,  Ey=0,  H=0,
                       omega=0, phi=0, T=0,
                       n=0, dt=0, alltime=0):
        self.program = "./flat"
        if sys.platform.startswith("win32"):
            self.program = "flat.exe"
        if dumping:
            self.program += " -d"
        self.args = "%s %e %e %e %e %e %e %e %e %e %d %e %e\n" % (m, Exc, Eyc, Hc, Ex, Ey, H, omega, phi, T, n, dt, alltime)

    def _parse(self, name, text):
        s = re.search(name + r"\s*=\s*(\S*)\s*\+/-\s*(\S*)", text)
        if s:
            mean, std = s.groups()
            return [float(mean), float(std)]
        return None

    def run(self, verbose=False):
        out, err = Popen([self.program], shell=True, stdin=PIPE, stdout=PIPE)\
                   .communicate(input=self.args.encode("utf8"))
        data = out.decode("utf8")
        if verbose:
            print(data)
        vx = self._parse("v_x", data)
        vy = self._parse("v_y", data)
        power = self._parse("power", data)
        tau = self._parse("tau", data)
        return {
            "v_x": vx,
            "v_y": vy,
            "power": power,
            "tau": tau
        }
