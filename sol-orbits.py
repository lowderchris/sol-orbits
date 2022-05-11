import astrospice
from sunpy.coordinates import HeliographicCarrington

from datetime import datetime, timedelta
import astropy.units as u
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from astropy.visualization import quantity_support
#quantity_support()

# Load the spice kernel
k_psp = astrospice.registry.get_kernels('psp', 'predict')
k_solo = astrospice.registry.get_kernels('solar orbiter', 'predict')
k_sta = astrospice.registry.get_kernels('stereo-a', 'predict')
#k_stb = astrospice.registry.get_kernels('stereo-b', 'recon')

# Generate some times
starttime = datetime(2022, 1, 1)
istarttime = starttime
endtime = starttime + timedelta(days=27.2753)
times = []
while istarttime < endtime:
    times.append(istarttime)
    istarttime += timedelta(hours=12)

# Generate some positions
c_psp = astrospice.generate_coords('SOLAR PROBE PLUS', times)
c_solo = astrospice.generate_coords('Solar orbiter', times)
c_sta = astrospice.generate_coords('Stereo ahead', times)
#c_stb = astrospice.generate_coords('Stereo behind', times)

# Convert coordinate systems
#ref_frame = HeliographicCarrington(observer='self')
#c_psp = c_psp.transform_to(ref_frame)
#c_solo = c_psp.transform_to(ref_frame)
#c_sta = c_psp.transform_to(ref_frame)
#c_stb = c_psp.transform_to(ref_frame)

# Convert to cartesian
x_psp = c_psp.cartesian
x_solo = c_solo.cartesian
x_sta = c_sta.cartesian
#x_stb = c_psp.cartesian

# Plot the orbit in 3d
times_float = [(t - times[0]).total_seconds() for t in times]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
kwargs = {'s': 3, 'c': times_float}
ax.scatter(x_psp.x.to(u.au), x_psp.y.to(u.au), x_psp.z.to(u.au), **kwargs)
ax.text(x_psp.x.to(u.au).value[-1], x_psp.y.to(u.au).value[-1], x_psp.z.to(u.au).value[-1], 'PSP')
ax.scatter(x_solo.x.to(u.au), x_solo.y.to(u.au), x_solo.z.to(u.au), **kwargs)
ax.text(x_solo.x.to(u.au).value[-1], x_solo.y.to(u.au).value[-1], x_solo.z.to(u.au).value[-1], 'SolO')
ax.scatter(x_sta.x.to(u.au), x_sta.y.to(u.au), x_sta.z.to(u.au), **kwargs)
ax.text(x_sta.x.to(u.au).value[-1], x_sta.y.to(u.au).value[-1], x_sta.z.to(u.au).value[-1], 'STA')
#ax.scatter(earth.x.to(u.au), earth.y.to(u.au), earth.z.to(u.au), **kwargs)
#ax.text(earth.x.to(u.au).value[-1], earth.y.to(u.au).value[-1], earth.z.to(u.au).value[-1], 'Earth')
ax.scatter(0,0,0, c='k')
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
ax.set_title(starttime.strftime("%Y-%m-%d %H:%M:%S") + ' --- ' + endtime.strftime("%Y-%m-%d %H:%M:%S"))

ax.plot([x_psp.x.to(u.au).value[0],0],[x_psp.y.to(u.au).value[0],0],[x_psp.z.to(u.au).value[0],0],'k', alpha=0.2)
ax.plot([x_psp.x.to(u.au).value[-1],0],[x_psp.y.to(u.au).value[-1],0],[x_psp.z.to(u.au).value[-1],0],'k', alpha=0.2)
ax.plot([x_solo.x.to(u.au).value[0],0],[x_solo.y.to(u.au).value[0],0],[x_solo.z.to(u.au).value[0],0],'k', alpha=0.2)
ax.plot([x_solo.x.to(u.au).value[-1],0],[x_solo.y.to(u.au).value[-1],0],[x_solo.z.to(u.au).value[-1],0],'k', alpha=0.2)
ax.plot([x_sta.x.to(u.au).value[0],0],[x_sta.y.to(u.au).value[0],0],[x_sta.z.to(u.au).value[0],0],'k', alpha=0.2)
ax.plot([x_sta.x.to(u.au).value[-1],0],[x_sta.y.to(u.au).value[-1],0],[x_sta.z.to(u.au).value[-1],0],'k', alpha=0.2)
#ax.plot([earth.x.to(u.au).value[0],0],[earth.y.to(u.au).value[0],0],[earth.z.to(u.au).value[0],0],'k', alpha=0.2)
#ax.plot([earth.x.to(u.au).value[-1],0],[earth.y.to(u.au).value[-1],0],[earth.z.to(u.au).value[-1],0],'k', alpha=0.2)

ax.set_xticks(ax.get_xticks()[::2])
ax.set_yticks(ax.get_yticks()[::2])
ax.set_zticks(ax.get_zticks()[::2])

tight_layout()

savefig('sol-orbits.pdf')
savefig('sol-orbits.png')

# TODO - Sort out the plot below once coordinate transformations have been completed

# Plot some orbit quantities
#fig, axs = plt.subplots(3, 1, sharex=True)
##axs[0].plot(earth.times, earth.r, 'k', label='Earth')
#axs[0].plot(times, c_psp.r, label='PSP')
#axs[0].plot(solo.times, solo.r, label='SolO')
#axs[0].plot(sta.times, sta.r, label='STA')
#axs[0].set_ylim(0, 1.1)
#axs[0].set_ylabel('r (AU)')
#axs[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#           ncol=4, mode="expand", borderaxespad=0.)

##axs[1].plot(earth.times, np.rad2deg(np.arcsin(earth.z / earth.r)), 'k')
#axs[1].plot(psp.times, np.rad2deg(np.arcsin(psp.z / psp.r)))
#axs[1].plot(solo.times, np.rad2deg(np.arcsin(solo.z / solo.r)))
#axs[1].plot(sta.times, np.rad2deg(np.arcsin(sta.z / sta.r)))
#axs[1].set_ylabel('Elevation (deg)')

##axs[2].plot(earth.times, earth.speed, 'k')
#axs[2].plot(psp.times, psp.speed)
#axs[2].plot(solo.times, solo.speed)
#axs[2].plot(sta.times, sta.speed)
#axs[2].set_ylabel('Speed (km/s)')

#fig.autofmt_xdate()

#tight_layout()

#savefig('sol-orbits-param.pdf')
