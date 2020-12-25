import numpy as np

y = np.array([2, 1, 5, 2])
x = np.array([1, 2, 3, 4])
m = np.ma.masked_where(y > 2, y)

print(m)
print(list(m))
print(np.ma.compressed(m))

print(y > 2)
print(np.where(y > 2))
