import numpy as np

test = np.load("test.npy")
print test.shape
print test[0, 0, 0]
for s in [0, 1, 2, 3]:
    print test[:, :, s]
