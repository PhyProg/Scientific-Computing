import numpy as np

x = np.array([[0],[1]])
y = np.array([[0],[3]])
z = []

y = np.c_[y, x]
print(y)
