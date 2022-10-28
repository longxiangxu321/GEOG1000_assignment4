import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

dev_x = [5000, 500000, 5000000, 50000000]

t_python = [0.038, 3.454, 34.946, 359.025]

t_cpp_release = [0.003, 0.29, 2.869, 29.239]

t_cpp_debug = [0.003, 0.029, 0.268, 2.708]

plt.plot(dev_x, t_python, label="Python")
plt.plot(dev_x, t_cpp_debug, label="C++ (Debug)")
plt.plot(dev_x, t_cpp_release, label="C++ (Release)")

plt.xlabel("Iteration")
plt.ylabel("Runtime")
plt.title("Runtime of nbody programs")

plt.legend()
plt.show()
