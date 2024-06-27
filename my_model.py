import math
import sys

import dnest4
import numpy as np
from scipy.special import gammaln

from data import Data
from lookup import Lookup
from my_distribution import MyDistribution


class MyModel:
    def __init__(
        self,
        num_dimensions=5,
        max_num_components=10,
        fixed=False,
        conditional_prior=None,
    ):
        self.num_dimensions = num_dimensions
        self.max_num_components = max_num_components
        self.fixed = fixed
        self.conditional_prior = conditional_prior or MyDistribution()
        self.components = []
        self.u_components = []

    def from_prior(self, rng=None):
        rng = rng or np.random.default_rng()
        self.components.clear()
        self.u_components.clear()

        num = (
            self.max_num_components
            if self.fixed
            else rng.integers(0, self.max_num_components + 1)
        )

        for _ in range(num):
            component = rng.random(self.num_dimensions)
            self.u_components.append(component.copy())
            self.conditional_prior.from_uniform(component)
            self.components.append(component)

        data_instance = Data.get_instance()
        data_instance.load("data/nu_oph.txt")
        background = (
            data_instance.get_y_min()
            + (data_instance.get_y_max() - data_instance.get_y_min()) * rng.random()
        )
        extra_sigma = math.exp(math.tan(math.pi * (0.97 * rng.random() - 0.485)))
        nu = math.exp(math.log(0.1) + math.log(1000.0) * rng.random())

        obj = self.serialize_objects()
        params = np.concatenate([obj, [background, extra_sigma, nu]])
        mu = self.calculate_mu(params)

        return np.concatenate([params, mu])

    def perturb(self, params, rng=None):
        logH = 0.0
        rng = rng or np.random.default_rng()
        if rng.random() <= 0.75:
            logH += self._perturb_components(rng)
            self.calculate_mu(params)
        elif rng.random() <= 0.5:
            params[-2] = np.log(params[-2])
            params[-2] = (np.arctan(params[-2]) / np.pi + 0.485) / 0.97
            params[-2] += rng.random()
            params[-2] = dnest4.wrap(params[-2], 0.0, 1.0)
            params[-2] = np.tan(np.pi * (0.97 * params[-2] - 0.485))
            params[-2] = np.exp(params[-2])

            params[-1] = np.log(params[-1])
            params[-1] += np.log(1000.0) * rng.random()
            params[-1] = dnest4.wrap(params[-1], np.log(0.1), np.log(1000.0))
            params[-1] = np.exp(params[-1])
        else:
            params[-3] += (
                Data.get_instance().get_y_max() - Data.get_instance().get_y_min()
            ) * rng.random()
            params[-3] = dnest4.wrap(
                params[-3],
                Data.get_instance().get_y_min(),
                Data.get_instance().get_y_max(),
            )
            self.calculate_mu(params)

        return logH

    def log_likelihood(self, params):
        y = Data.get_instance().get_y()
        sig = Data.get_instance().get_sig()
        extra_sigma = params[-2]
        nu = params[-1]

        logL = 0.0
        for i in range(len(y)):
            var = sig[i] ** 2 + extra_sigma**2
            logL += gammaln(0.5 * (nu + 1.0)) - gammaln(0.5 * nu)
            logL += -0.5 * np.log(np.pi * nu) - 0.5 * np.log(var)
            logL += (
                -0.5 * (nu + 1.0) * np.log(1.0 + (y[i] - params[-1]) ** 2 / (var * nu))
            )

        return logL

    def print(self, params, out=sys.stdout):
        t = Data.get_instance().get_t()
        t_min, t_max = min(t), max(t)
        finer_t = np.linspace(t_min, t_max, 1000)

        components = self.deserialize_objects(params[:-3])
        signal = np.full_like(finer_t, params[-3])

        for component in components:
            T = math.exp(component[0])
            A = component[1]
            phi = component[2]
            try:
                v0 = math.sqrt(1.0 - component[3])
            except:
                v0 = 0
            viewing_angle = component[4]

            arg = 2.0 * math.pi * finer_t / T + phi
            evaluations = Lookup.get_instance().evaluate(arg, v0, viewing_angle)
            signal += A * np.array(evaluations)

        np.savetxt(out, signal)
        out.write(f"{params[-2]} {params[-1]} ")
        self._print_components(out)
        out.write("\n")

    def calculate_mu(self, params):
        t = np.array(Data.get_instance().get_t())
        components = self.deserialize_objects(params[:-3])
        background = params[-3]

        mu = np.full(len(t), background)
        for component in components:
            T = math.exp(component[0])
            A = component[1]
            phi = component[2]
            try:
                v0 = math.sqrt(1.0 - component[3])
            except:
                v0 = 0
            viewing_angle = component[4]

            arg = 2.0 * math.pi * t / T + phi
            evaluations = Lookup.get_instance().evaluate(arg, v0, viewing_angle)
            mu += A * np.array(evaluations)
        return mu

    def serialize_objects(self):
          return np.hstack([np.array(component) for component in self.components]).flatten()

    def deserialize_objects(self, array):
        components = []
        num_components = len(array) // self.num_dimensions
        for i in range(num_components):
            components.append(array[i * self.num_dimensions : (i + 1) * self.num_dimensions])
        return components

    def _perturb_components(self, rng):
        logH = 0.0

        if rng.random() < 0.5 and not self.fixed:
            if rng.random() < 0.5 and len(self.components) < self.max_num_components:
                component = rng.random(self.num_dimensions)
                self.u_components.append(component.copy())
                self.conditional_prior.from_uniform(component)
                self.components.append(component)
            elif len(self.components) > 0:
                idx = rng.integers(len(self.components))
                self.components.pop(idx)
                self.u_components.pop(idx)
            logH = 0.0
        else:
            if len(self.components) > 0:
                idx = rng.integers(len(self.components))
                self.u_components[idx] += rng.random(self.num_dimensions)
                self.u_components[idx] %= 1.0
                self.components[idx] = self.u_components[idx].copy()
                self.conditional_prior.from_uniform(self.components[idx])
            logH = 0.0

        return logH

    def _print_components(self, out=sys.stdout):
        for component in self.components:
            out.write(" ".join(map(str, component)) + " ")
        out.write("\n")

    def description(self):
        return ""
