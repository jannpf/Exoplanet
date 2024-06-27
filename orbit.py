import math
import numpy as np
import sys


class Orbit:
    def __init__(self):
        self.vx = []
        self.vy = []

    def load(self, filename):
        try:
            with open(filename, "r") as fin:
                lines = fin.readlines()
                for line in lines:
                    temp1, temp2 = map(float, line.split())
                    self.vx.append(temp1)
                    self.vy.append(temp2)
            print(f"# Loaded {len(self.vx)} points from orbit file {filename}.")
        except IOError:
            print(f"# ERROR: Couldn't open file {filename}.")
            sys.exit(1)

    def evaluate(self, arg_to_cos, viewing_angle):
        C = math.cos(viewing_angle)
        S = math.sin(viewing_angle)
        vr = [C * x + S * y for x, y in zip(self.vx, self.vy)]

        vrmax = max(vr)
        vrmin = min(vr)
        highest = vr.index(vrmax)

        # Normalize vr
        vr = [2 * (v - vrmin) / (vrmax - vrmin) - 1 for v in vr]

        result = list(arg_to_cos)
        cc = 2.0 * math.pi * highest / len(vr)
        for i in range(len(result)):
            index = int(
                len(vr)
                * (self._mod(arg_to_cos[i] + cc, 2.0 * math.pi) / (2.0 * math.pi))
            )
            result[i] = vr[index]

        return result

    def _mod(self, a, b):
        return a - b * math.floor(a / b)

    @staticmethod
    def test():
        o = Orbit()
        o.load("artifacts/orbits0.710.dat")

        t = [tt for tt in np.arange(-10, 10.01, 0.01)]
        y = o.evaluate(t, 1.0)
        for ti, yi in zip(t, y):
            print(ti, yi)


if __name__ == "__main__":
    Orbit.test()
