"""
The script simulates the orbital motion of a body under a central gravitational
force for a range of initial velocities. It calculates the position and velocity
of the body at each timestep, saving the velocities to files named based on the
initial velocity. These files are later used to analyze the orbital characteristics and radial velocities.

Some details:
1. Iterates over 121 values of initial velocity (v0) ranging from 0.4 to 1.0.
2. Calculates velocities for each timestep for each of the initial velocities.
3. Saves data for each initial velocity to a separate file.

This data is later used in orbit_lookup.py.
"""

import pylab


def main():
    for v0 in pylab.linspace(0.4, 1.0, 121):
        # Initial position
        pos = pylab.array([1.0, 0.0])

        # Initial velocity
        vel = pylab.array([0.0, v0])

        # Timestep
        dt = 0.0005
        tfinal = 10.0
        steps = int(tfinal / dt)

        keep = pylab.empty((steps, 2))

        pylab.ion()
        for i in range(0, steps):
            phi1 = pylab.arctan2(pos[1], pos[0])

            pos += 0.5 * dt * vel
            accel = -pos / (pos[0] ** 2 + pos[1] ** 2) ** 1.5
            vel += dt * accel
            pos += 0.5 * dt * vel
            keep[i, :] = vel

            phi2 = pylab.arctan2(pos[1], pos[0])

            # 	  if i%50 == 0:
            # 	    hold(False)
            # 	    plot(keep[0:(i+1), 0], keep[0:(i+1), 1], 'bo')
            # 	    hold(True)
            # 	    plot(0, 0, 'r*')
            # 	    axis('equal')
            # 	    draw()

            if phi1 < 0.0 and phi2 > 0.0:
                keep = keep[0 : (i + 1), :]
                break

        pylab.savetxt("../artifacts/" + "orbits{v:.3f}".format(v=v0) + ".dat", keep)

    pylab.ioff()
    pylab.show()


if __name__ == "__main__":
    main()
