import numpy as np
import matplotlib.pyplot as plt
#
a = 1.4
b = 0.3
X = [(0, 0)]
for n in range(100000):
    X.append((X[-1][1] + 1 - a * X[-1][0]**2, b * X[-1][0]))
#
#fig = plt.figure()
#plt.plot([x[0] for x in X],[x[1] for x in X],'.b')

nbins = 200
H, xedges, yedges = np.histogram2d(
    [x[0] for x in X], [x[1] for x in X], bins=nbins)

# H needs to be rotated and flipped
H = np.rot90(H)
H = np.flipud(H)

# Mask zeros
Hmasked = np.ma.masked_where(H == 0, H)  # Mask pixels with a value of zero

# Plot 2D histogram using pcolor
fig2 = plt.figure()
plt.pcolormesh(xedges, yedges, Hmasked)
plt.xlabel('x')
plt.ylabel('y')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Counts')
plt.show()