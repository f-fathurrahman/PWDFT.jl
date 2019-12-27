import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt("DATA_c2c.txt")

plt.clf()
plt.plot(dat[:,0]/1e6, dat[:,1], marker="o", label="CPU timing")
plt.plot(dat[:,0]/1e6, dat[:,2], marker="o", label="GPU timing")
plt.ylabel("Timing (ms)")
plt.xlabel("Data size (x$10^6$)")
plt.grid()
plt.legend()
plt.title("Complex to complex 3D FFT operation")
plt.savefig("DATA_c2c.pdf")

timingCPU = dat[:,1]
timingGPU = dat[:,2]
for i in range(len(timingCPU)):
    print("%f" % (timingCPU[i]/timingGPU[i]))
