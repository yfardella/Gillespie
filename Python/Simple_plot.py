import numpy as np
import matplotlib.pyplot as plt

data_r = np.loadtxt('.\\c.txt')
data_time = np.loadtxt('.\\time.txt')

plt.grid(True)
plt.plot(data_time, data_r)
plt.plot(data_time, np.ones(data_time.shape) * 105.5/2.)

plt.ylim([0, 200])
plt.xlim([0, 10])
plt.xlabel(r'time')
plt.ylabel(r'receptors')
# plt.ylim([0,5])
plt.show()