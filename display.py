from pylab import *
import os
import dnest4.classic as dn4

rc("font", size=14, family="serif", serif="Computer Sans")
rc("text", usetex=True)

star = 'hd191939'
# star = 'hip14810'
# star = 'corot7'

data = loadtxt(f'{star}.txt')
# truth = loadtxt('fake_data_like_nuoph.truth')
posterior_sample = atleast_2d(dn4.my_loadtxt('posterior_sample.txt'))

width=0.5
hist(posterior_sample[:,1007], bins=arange(0, 11)-0.5*width, width=width, color="k", alpha=0.2)
xlabel('Number of Planets')
ylabel('Number of Posterior Samples')
xlim([-0.5, 10.5])
get_current_fig_manager().set_window_title(os.path.dirname(os.path.realpath(__file__)))
show()

T = posterior_sample[:,1008:1018]
A = posterior_sample[:,1018:1028]
E = posterior_sample[:,1038:1048]
which = T != 0
T = T[which].flatten()
A = A[which].flatten()
E = E[which].flatten()
# Trim
#s = sort(T)
#left, middle, right = s[0.25*len(s)], s[0.5*len(s)], s[0.75*len(s)]
#iqr = right - left
#s = s[logical_and(s > middle - 5*iqr, s < middle + 5*iqr)]

hist(T/log(10.), 500, alpha=0.4, color="k")
xlabel(r'$\log_{10}$(Period/days)')
xlim([0, 5])

true_periods = {
  'hd191939': [
    8.8803256, # hd191939
    28.579743, # hd191939
    38.353037, # hd191939
    101.12, # hd191939
    284, # hd191939
    2200 # hd191939
  ],
  'hip14810': [
    6.67, # hip 14810
    147.7, # hip 14810
    952, # hip 14810
  ],
  'corot7': [
    0.853592, # corot-7
    3.697, # corot-7
    8.966 # corot-7
  ]
}

true_eccs = {
  'hd191939': [
    0.031, # hd191939
    0.034, # hd191939
    0.031, # hd191939
    0.031, # hd191939
    0.030, # hd191939
    nan # 2200 # hd191939 # no data
  ],
  'hip14810': [
    0.14399, # hip 14810
    0.1566, # hip 14810
    0.185, # hip 14810
  ],
  'corot7': [
    0, # corot-7
    0, # corot-7
    0, # corot-7
  ]
}

if star in true_periods.keys():
  for p in true_periods[star]:
    axvline(log(p)/log(10.), color='g')
ylabel('Number of Posterior Samples')
get_current_fig_manager().set_window_title(os.path.dirname(os.path.realpath(__file__)))
show()

subplot(2,1,1)
# plot(truth[1008:1008 + int(truth[1007])]/log(10.), log10(truth[1018:1018 + int(truth[1007])]), 'ko', markersize=7, alpha=0.5)
xlim([0, 5])
ylim([-1, 4])
ylabel(r'$\log_{10}$[Amplitude (m/s)$]$')
plot(T/log(10.), log10(A), 'g.', markersize=1)

subplot(2,1,2)
plot(log(true_periods[star])/log(10.), true_eccs[star], 'ko', markersize=7, alpha=0.5)
xlim([0, 5])
xlabel(r'$\log_{10}$(Period/days)')
ylabel('Eccentricity')
plot(T/log(10.), E, 'g.', markersize=1)
get_current_fig_manager().set_window_title(os.path.dirname(os.path.realpath(__file__)))
show()

data[:,0] -= data[:,0].min()
t = linspace(data[:,0].min(), data[:,0].max(), 1000)

saveFrames = False # For making movies
if saveFrames:
  os.system('rm Frames/*.png')

for i in range(0, posterior_sample.shape[0]):
  clf()
  errorbar(data[:,0], data[:,1], fmt='k.', yerr=data[:,2])
  plot(t, posterior_sample[i, 0:1000], 'g')
  xlim([-0.05*data[:,0].max(), 1.05*data[:,0].max()])
  ylim([-1.5*max(abs(data[:,1])), 1.5*max(abs(data[:,1]))])
  #axhline(0., color='k')
  xlabel('Time (days)', fontsize=16)
  ylabel('Radial Velocity (m/s)', fontsize=16)
  if saveFrames:
    savefig('Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
    print('Frames/' + '%0.4d'%(i+1) + '.png')


get_current_fig_manager().set_window_title(os.path.dirname(os.path.realpath(__file__)))
show()
