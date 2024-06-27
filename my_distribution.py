import math
import sys
import dnest4


class MyDistribution:
    def __init__(self):
        self.center = 0.0
        self.width = 0.0
        self.mu = 0.0

    def from_prior(self, r):
        # Cauchy prior centered on 5.901 = log(365 days).
        self.center = 5.901 + math.tan(math.pi * (0.97 * dnest4.rand() - 0.485))
        self.width = 0.1 + 2.9 * dnest4.rand()
        self.mu = math.exp(math.tan(math.pi * (0.97 * dnest4.rand() - 0.485)))

    def perturb_hyperparameters(self, rng):
        logH = 0.0

        which = int(dnest4.rand(3))

        if which == 0:
            self.center = (math.atan(self.center - 5.901) / math.pi + 0.485) / 0.97
            self.center += dnest4.randh()
            self.center = 5.901 + math.tan(math.pi * (0.97 * self.center - 0.485))
        elif which == 1:
            self.width += 2.9 * dnest4.randh()
            self.width = self._wrap(self.width, 0.1, 3.0)
        else:
            self.mu = math.log(self.mu)
            self.mu = (math.atan(self.mu) / math.pi + 0.485) / 0.97
            self.mu += dnest4.randh()
            self.mu = math.tan(math.pi * (0.97 * self.mu - 0.485))
            self.mu = math.exp(self.mu)

        return logH

    def log_pdf(self, vec):
        # vec[0] = "position" (log-period)
        # vec[1] = amplitude
        # vec[2] = phase
        # vec[3] = v0
        # vec[4] = viewing angle
        if (
            vec[1] < 0.0
            or vec[2] < 0.0
            or vec[2] > 2.0 * math.pi
            or vec[3] < 0.0
            or vec[3] > 0.8189776
            or vec[4] < 0.0
            or vec[4] > 2.0 * math.pi
        ):
            return -1e300

        return (
            -math.log(2.0 * self.width)
            - abs(vec[0] - self.center) / self.width
            - math.log(self.mu)
            - vec[1] / self.mu
            + 2.1 * math.log(1.0 - vec[3] / 0.995)
        )

    def from_uniform(self, vec):
        if vec[0] < 0.5:
            self.center = self.center + self.width * math.log(2.0 * vec[0])
        else:
            self.center = self.center - self.width * math.log(2.0 - 2.0 * vec[0])
        vec[1] = -self.mu * math.log(1.0 - vec[1])
        vec[2] = 2.0 * math.pi * vec[2]
        vec[3] = 1.0 - math.pow(1.0 - 0.995 * vec[3], 1.0 / 3.1)
        vec[4] = 2.0 * math.pi * vec[4]

    def to_uniform(self, vec):
        if vec[0] < self.center:
            vec[0] = 0.5 * math.exp((vec[0] - self.center) / self.width)
        else:
            vec[0] = 1.0 - 0.5 * math.exp((self.center - vec[0]) / self.width)
        vec[1] = 1.0 - math.exp(-vec[1] / self.mu)
        vec[2] = vec[2] / (2.0 * math.pi)
        vec[3] = 1.0 - math.pow(1.0 - vec[3] / 0.995, 3.1)
        vec[4] = vec[4] / (2.0 * math.pi)

    def print(self, out=sys.stdout):
        print(self.center, self.width, self.mu, file=out)

    def _wrap(self, x, xmin, xmax):
        width = xmax - xmin
        while x < xmin:
            x += width
        while x >= xmax:
            x -= width
        return x


if __name__ == "__main__":
    import random

    # dnest4 = random.Random()
    distribution = MyDistribution()

    distribution.from_prior(None)
    print("Parameters after initialization:")
    distribution.print()

    vec = [0.3, 0.5, 1.0, 0.2, 1.5]
    print(f"Log PDF for vector {vec}: {distribution.log_pdf(vec)}")
