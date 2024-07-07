import os

import dnest4.classic as dn4
import pylab

pylab.rc("font", size=14, family="serif", serif="Computer Sans")
pylab.rc("text")

star = "nu_oph"
true_periods = {
    "nu_oph": [530.73, 3188.95],  # b, c
}

true_eccs = {
    "nu_oph": [0.124, 0.178],  # b, c
}

data = pylab.loadtxt("../../nu_oph.txt")
posterior_sample = pylab.atleast_2d(dn4.my_loadtxt("posterior_sample.txt"))

width = 0.5
pylab.hist(
    posterior_sample[:, 1007],
    bins=pylab.arange(0, 12) - 0.5 * width,
    width=width,
    color="k",
    alpha=0.2,
)
pylab.xlabel("Number of Planets")
pylab.ylabel("Number of Posterior Samples")
pylab.xlim([-0.5, 10.5])
pylab.show()

T = posterior_sample[:, 1008:1018]
A = posterior_sample[:, 1018:1028]
E = posterior_sample[:, 1038:1048]

which = T != 0
T = T[which].flatten()
A = A[which].flatten()
E = E[which].flatten()
# Trim
# s = sort(T)
# left, middle, right = s[0.25*len(s)], s[0.5*len(s)], s[0.75*len(s)]
# iqr = right - left
# s = s[logical_and(s > middle - 5*iqr, s < middle + 5*iqr)]

pylab.hist(T / pylab.log(10.0), bins=500, alpha=0.2, color="k")
pylab.xlabel(r"$\log_{10}$(Period/days)")
pylab.xlim([1, 5])
pylab.ylabel("Number of Posterior Samples")
if star in true_periods.keys():
    for p in true_periods[star]:
        pylab.axvline(pylab.log(p) / pylab.log(10.0), color="g", alpha=0.5)
pylab.show()

pylab.subplot(2, 1, 1)
pylab.xlim([1, 4])
pylab.ylim([-1, 3])
pylab.ylabel(r"$\log_{10}$[Amplitude (m/s)$]$")
pylab.plot(T / pylab.log(10.0), pylab.log10(A), "g.", markersize=1)

pylab.subplot(2, 1, 2)
pylab.plot(
    pylab.log(true_periods[star]) / pylab.log(10.0), true_eccs[star], "ko", markersize=7, color='r'
)
pylab.xlim([1, 5])
pylab.xlabel(r"$\log_{10}$(Period/days)")
pylab.ylabel("Eccentricity")
pylab.plot(T / pylab.log(10.0), E, "g.", markersize=1)
pylab.show()

data[:, 0] -= data[:, 0].min()
t = pylab.linspace(data[:, 0].min(), data[:, 0].max(), 1000)

saveFrames = False  # For making movies
if saveFrames:
    os.system("rm Frames/*.png")

for i in range(0, posterior_sample.shape[0]):
    pylab.clf()
    pylab.errorbar(data[:, 0], data[:, 1], fmt="k.", yerr=data[:, 2])
    pylab.plot(t, posterior_sample[i, 0:1000], "g")
    pylab.xlim([-0.05 * data[:, 0].max(), 1.05 * data[:, 0].max()])
    pylab.ylim([-1.5 * max(abs(data[:, 1])), 1.5 * max(abs(data[:, 1]))])
    # axhline(0., color='k')
    pylab.xlabel("Time (days)", fontsize=16)
    pylab.ylabel("Radial Velocity (m/s)", fontsize=16)
    if saveFrames:
        pylab.savefig("Frames/" + "%0.4d" % (i + 1) + ".png", bbox_inches="tight")
        print("Frames/" + "%0.4d" % (i + 1) + ".png")

pylab.show()

