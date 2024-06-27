"""
This combined script facilitates the visualization and comparison
of simulated orbital data with the output from the kepler.py script,
likely to study how different orbital parameters and viewing
angles affect the observed radial velocity.


What it does:
1. Calls kepler.py to simulate "true" keplerian orbit
(../artifacts/kepler.txt, containing time-rv pairs).
2. Loads trajectories created by orbits.py. These trajectories
contain time-velocity data for different initial velocities.
Velocities are estimated based on system parameters such as
eccentricity, etc.
3. Runs evaluate.py which evaluates radial velocity based
on trajectory data.
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # Execute the kepler.py script using system command
    subprocess.run(["python3", "kepler.py"])

    # Load all the files
    trajectories = []
    for i in range(1, 122):
        v = 0.4 + 0.005 * (i - 1)
        print(v)
        filename = f"../artifacts/orbits{v:.3f}.dat"
        try:
            trajectory = np.loadtxt(filename)
        except FileNotFoundError:
            subprocess.run(["python3", "orbits.py"])
            trajectory = np.loadtxt(filename)
        trajectories.append(trajectory)

    def evaluate(arg_to_sin, v=0.75, viewing_angle=0.):
        """
        Which trajectory to use?
        Computes the radial velocity based on the selected trajectory, velocity v, and viewing angle.
        Determines the appropriate trajectory to use based on v.
        Calculates the radial velocity by combining the trajectory data and normalizes it.
        Finds the index corresponding to the maximum radial velocity and adjusts the phase to compute the final index.
        Returns the evaluated radial velocity value.
        """

        orbit = int(np.floor((v - 0.4) / 0.005)) + 1
        N = trajectories[orbit - 1].shape[0]

        radial_velocity = np.cos(viewing_angle) * trajectories[orbit - 1][:, 0] + np.sin(viewing_angle) * trajectories[orbit - 1][:, 1]
        radial_velocity /= np.max(np.abs(radial_velocity))
        radial_velocity = (radial_velocity - np.min(radial_velocity)) / (np.max(radial_velocity) - np.min(radial_velocity))
        radial_velocity = 2 * radial_velocity - 1
        closest_to_max = np.where(radial_velocity == np.max(radial_velocity))[0][0]

        cc = 2 * np.pi * (closest_to_max - 1) / N
        index = 1 + N * ((arg_to_sin + cc) % (2 * np.pi)) / (2 * np.pi)

        return radial_velocity[int(index)]

    t = np.arange(0, 10, 0.01)
    y = [evaluate(2 * np.pi * ti / 5 + 2.2, v=np.sqrt(1. - 0.5), viewing_angle=-2.7) for ti in t]
    plt.plot(t, y)

    data = pd.read_csv("../artifacts/kepler.txt", sep=" ", header=None)
    plt.plot(data.iloc[:, 0], data.iloc[:, 1])

    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.title('Plot of y vs t and data from kepler.txt')
    plt.legend(['y vs t', 'Data from kepler.txt'])
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
