import numpy as np

x = np.array([0.2, 6.4, 3.0, 1.6, 10.0])
bins = np.array([0.0, 1.0, 2.5, 4.0, 10.0])

# bins[i-1] <= x < bins[i]
# lower: [0, 3, 2, 1, 4] (i-1)
# upper: [1, 4, 3, 2, 5] (i)
# last index is out of bin range
inds = np.digitize(x, bins)
print(inds)

# bins[i-1] < x <= bins[i]
# lower: [0, 3, 2, 1, 3] (i-1)
# upper: [1, 4, 3, 2, 4] (i)
inds = np.digitize(x, bins, right=True)
print(inds)
