from viewRSoftData import *
import matplotlib.pyplot as plt
import numpy as np
from lanternfiber import lanternfiber
plt.ion()

outpowers = np.array([0.11047041, 0.14290185, 0.32011984, 0.19727427, 0.18515128,
              0.15496   , 0.30159726, 0.25067494, 0.1288779 , 0.10380741,
              0.43701958, 0.09005142, 0.25019437, 0.19806124, 0.22968186,
              0.29892668, 0.19674209, 0.03508904, 0.21348379])

predpowers = np.array([0.10773016, 0.30387031, 0.13247846, 0.31291402, 0.13743055,
              0.18379147, 0.17573258, 0.1069473 , 0.44455322, 0.10866876,
              0.13136371, 0.27050398, 0.22000963, 0.04617284, 0.20052526,
              0.28847195, 0.23507836, 0.20093922, 0.25685981])


outpow_sort = np.argsort(outpowers)
predpow_sort = np.argsort(predpowers)

new_pred_inds = np.zeros_like(predpow_sort)
for k in range(len(outpow_sort)):
    # new_pred_inds[k] = np.where(predpow_sort == outpow_sort[k])[0]
    new_pred_inds[outpow_sort[k]] = predpow_sort[k]


plt.plot(outpowers)
plt.plot(predpowers[new_pred_inds])

