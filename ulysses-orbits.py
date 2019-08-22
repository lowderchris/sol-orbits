import heliopy.data.spice as spicedata
import heliopy.spice as spice
from datetime import datetime, timedelta
import astropy.units as u
import numpy as np

# Load the spice kernel
kernels = spicedata.get_kernel('ulysses')
spice.furnish(kernels)
uly = spice.Trajectory('ulysses')

# Generate some times
# Add support for specification of a CR
starttime = datetime(1990, 10, 7)
endtime = starttime + timedelta(days=6.2*365)
times = []
while starttime < endtime:
    times.append(starttime)
    starttime += timedelta(hours=24)

# Generate some positions
uly.generate_positions(times, 'Sun', 'ECLIPJ2000')
uly.change_units(u.au)

# Plot the orbit in 3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.visualization import quantity_support
quantity_support()

# Generate a set of timestamps to color the orbits by
times_float = [(t - uly.times[0]).total_seconds() for t in uly.times]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
kwargs = {'s': 3, 'c': times_float}
ax.scatter(uly.x, uly.y, uly.z, **kwargs)
ax.scatter(0,0,0, c='k')
ax.set_xlim(-5.5, 5.5)
ax.set_ylim(-5.5, 5.5)
ax.set_zlim(-5.5, 5.5)
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
ax.set_title(starttime.strftime("%Y-%m-%d %H:%M:%S") + ' --- ' + endtime.strftime("%Y-%m-%d %H:%M:%S"))


# Plot some orbit quantities
elevation = np.rad2deg(np.arcsin(uly.z / uly.r))

fig, axs = plt.subplots(3, 1, sharex=True)
axs[0].plot(uly.times, uly.r)
axs[0].set_ylim(0, 5.5)
axs[0].set_ylabel('r (AU)')

axs[1].plot(uly.times, elevation)
axs[1].set_ylabel('Elevation (deg)')

axs[2].plot(uly.times, uly.speed)
axs[2].set_ylabel('Speed (km/s)')

plt.show()
