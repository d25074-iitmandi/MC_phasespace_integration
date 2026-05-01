import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../results/scaling.txt")
threads = data[:,0]
time = data[:,1]
speedup = data[:,2]

plt.plot(threads, speedup, marker='o')
plt.xlabel("Threads")
plt.ylabel("Speedup")
plt.title("OpenMP Scaling")
plt.grid()
plt.show()
