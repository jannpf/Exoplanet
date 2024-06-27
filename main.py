import subprocess

import dnest4

from data import Data
from lookup import Lookup
from my_model import MyModel


def main():
    # Load the data and lookup tables
    data = Data.get_instance().load("data/fake_data_like_nuoph.txt")
    lookup = Lookup.get_instance().load()

    # Start the DNest4 framework with MyModel
    model = MyModel()
    sampler = dnest4.DNest4Sampler(
        model, backend=dnest4.backends.CSVBackend(".", sep=" ")
    )
    gen = sampler.sample(
        max_num_levels=30,
        num_steps=1000,
        new_level_interval=50000,
        num_per_step=10000,
        thread_steps=100,
        num_particles=5,
        lam=10,
        beta=100,
        seed=1234,
    )

    # Do the sampling (one iteration here = one particle save)
    # breakpoint()
    for i, sample in enumerate(gen):
        print("# Saved {k} particles.".format(k=(i + 1)))

if __name__ == "__main__":
    main()
    subprocess.run("python3 ./helpers/showresults.py")
