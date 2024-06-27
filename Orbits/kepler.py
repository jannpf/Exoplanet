"""
This code essentially simulates the radial velocity curve
of an orbiting body, which is a fundamental technique in detecting
and studying exoplanets.


Details:
1. Set up time and orbital parameters.
2. Compute eccentric anomaly E and true anomaly f to obtain rv from E.
3. Plot and save the data to kepler.txt.

The generated file kepler.txt contains two columns: the first column is time values 𝑡
and the second column is the corresponding radial velocity 𝑣.
This data likely represents the radial velocity of a star or planet over time,
which is a common way to study exoplanetary systems using Kepler's laws of planetary motion.
"""


import pylab


def main():
    t = pylab.linspace(0, 10, 10001)

    e = 0.5  # Eccentricity
    P = 5.0  # Period
    K = 1.0  # Semi-amplitude
    phi = 1.8  # Longitude of ascending node
    wbar = 2.7  # Longitude of periastron
    omega = 2.0 * pylab.pi / P

    # Solve for E
    def update(E, t):
        E = E - (E - e * pylab.sin(E) - omega * t - phi) / (1.0 - e * pylab.cos(E))
        return E

    pylab.ion()
    E = omega * t + phi
    for i in range(0, 20):
        pylab.plot(t, E)
        E = update(E, t)
        pylab.title(f"{i + 1}")
        pylab.draw()

    pylab.ioff()
    pylab.show()

    E = pylab.mod(E, 2 * pylab.pi)
    cosf = (pylab.cos(E) - e) / (1.0 - e * pylab.cos(E))
    f = pylab.arccos(cosf)
    f[E > pylab.pi] = 2 * pylab.pi - f[E > pylab.pi]
    v = K * (pylab.sin(f + wbar) + e * pylab.sin(wbar))
    v -= (v.min() + v.max()) / 2
    print(v.max() - v.min())
    pylab.plot(t, v)

    data = pylab.empty((10001, 2))
    data[:, 0], data[:, 1] = t, v
    pylab.savetxt("../artifacts/kepler.txt", data)
    pylab.show()

    # EQUIVALENCE
    # kepler.py			Old code (but using cosine alignment)
    # (phi, wbar) = (0, pi/2)		(phase, viewing angle) = (0, -pi/2)
    # (phi, wbar) = (0, pi)			(phase, viewing angle) = (0.6, -pi)
    # (phi, wbar) = (0, 0.45)		(phase, viewing angle) = (5.9, -0.45)
    # (phi, wbar) = (1.1, 0.45)		(phase, viewing angle) = (0.7, -0.45)
    # (phi, wbar) = (1.1, 0.45)		(phase, viewing angle) = (0.7, -0.45)
    # (phi, wbar) = (1.8, 1.9)		(phase, viewing angle) = (1.9, -1.9)
    # (phi, wbar) = (1.8, 2.7)		(phase, viewing angle) = (2.2, -2.7)


if __name__ == "__main__":
    main()
