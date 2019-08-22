import heliopy.data.spice as spicedata
import heliopy.spice as spice
from datetime import datetime, timedelta
import astropy.units as u
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.visualization import quantity_support
quantity_support()

# Load the spice kernel
kernels = spicedata.get_kernel('psp')
kernels += spicedata.get_kernel('psp_pred')
kernels += spicedata.get_kernel('solo_2020')
kernels += spicedata.get_kernel('stereo_a')
kernels += spicedata.get_kernel('stereo_a_pred')
kernels += spicedata.get_kernel('planet_trajectories')
spice.furnish(kernels)
psp = spice.Trajectory('SPP')
solo = spice.Trajectory('solo')
sta = spice.Trajectory('stereo ahead')
earth = spice.Trajectory('earth')

# Generate some times
starttime = datetime(2020, 9, 15)
istarttime = starttime
endtime = starttime + timedelta(days=27.2753)
times = []
while istarttime < endtime:
    times.append(istarttime)
    istarttime += timedelta(hours=12)

# Generate some positions
psp.generate_positions(times, 'Sun', 'ECLIPJ2000')
psp.change_units(u.au)
solo.generate_positions(times, 'Sun', 'ECLIPJ2000')
solo.change_units(u.au)
sta.generate_positions(times, 'Sun', 'ECLIPJ2000')
sta.change_units(u.au)
earth.generate_positions(times, 'Sun', 'ECLIPJ2000')
earth.change_units(u.au)

# Plot the orbit in 3d
times_float = [(t - psp.times[0]).total_seconds() for t in psp.times]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
kwargs = {'s': 3, 'c': times_float}
ax.scatter(psp.x, psp.y, psp.z, **kwargs)
ax.text(psp.x.value[-1], psp.y.value[-1], psp.z.value[-1], 'PSP')
ax.scatter(solo.x, solo.y, solo.z, **kwargs)
ax.text(solo.x.value[-1], solo.y.value[-1], solo.z.value[-1], 'SolO')
ax.scatter(sta.x, sta.y, sta.z, **kwargs)
ax.text(sta.x.value[-1], sta.y.value[-1], sta.z.value[-1], 'STA')
ax.scatter(earth.x, earth.y, earth.z, **kwargs)
ax.text(earth.x.value[-1], earth.y.value[-1], earth.z.value[-1], 'Earth')
ax.scatter(0,0,0, c='k')
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_zlim(-1, 1)
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
ax.set_title(starttime.strftime("%Y-%m-%d %H:%M:%S") + ' --- ' + endtime.strftime("%Y-%m-%d %H:%M:%S"))

ax.plot([psp.x.value[0],0],[psp.y.value[0],0],[psp.z.value[0],0],'k', alpha=0.2)
ax.plot([psp.x.value[-1],0],[psp.y.value[-1],0],[psp.z.value[-1],0],'k', alpha=0.2)
ax.plot([solo.x.value[0],0],[solo.y.value[0],0],[solo.z.value[0],0],'k', alpha=0.2)
ax.plot([solo.x.value[-1],0],[solo.y.value[-1],0],[solo.z.value[-1],0],'k', alpha=0.2)
ax.plot([sta.x.value[0],0],[sta.y.value[0],0],[sta.z.value[0],0],'k', alpha=0.2)
ax.plot([sta.x.value[-1],0],[sta.y.value[-1],0],[sta.z.value[-1],0],'k', alpha=0.2)
ax.plot([earth.x.value[0],0],[earth.y.value[0],0],[earth.z.value[0],0],'k', alpha=0.2)
ax.plot([earth.x.value[-1],0],[earth.y.value[-1],0],[earth.z.value[-1],0],'k', alpha=0.2)

ax.set_xticks(ax.get_xticks()[::2])
ax.set_yticks(ax.get_yticks()[::2])
ax.set_zticks(ax.get_zticks()[::2])

tight_layout()

savefig('sol-orbits.pdf')
savefig('sol-orbits.png')

# Plot some orbit quantities
fig, axs = plt.subplots(3, 1, sharex=True)
axs[0].plot(earth.times, earth.r, 'k', label='Earth')
axs[0].plot(psp.times, psp.r, label='PSP')
axs[0].plot(solo.times, solo.r, label='SolO')
axs[0].plot(sta.times, sta.r, label='STA')
axs[0].set_ylim(0, 1.1)
axs[0].set_ylabel('r (AU)')
axs[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)

axs[1].plot(earth.times, np.rad2deg(np.arcsin(earth.z / earth.r)), 'k')
axs[1].plot(psp.times, np.rad2deg(np.arcsin(psp.z / psp.r)))
axs[1].plot(solo.times, np.rad2deg(np.arcsin(solo.z / solo.r)))
axs[1].plot(sta.times, np.rad2deg(np.arcsin(sta.z / sta.r)))
axs[1].set_ylabel('Elevation (deg)')

axs[2].plot(earth.times, earth.speed, 'k')
axs[2].plot(psp.times, psp.speed)
axs[2].plot(solo.times, solo.speed)
axs[2].plot(sta.times, sta.speed)
axs[2].set_ylabel('Speed (km/s)')

fig.autofmt_xdate()

tight_layout()

savefig('sol-orbits-param.pdf')
