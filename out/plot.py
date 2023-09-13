import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./P_rho_G_sv_1.txt")
data1 = np.loadtxt("./P_rho_G_sv_2.txt")
data2 = np.loadtxt("./P_rho_G_sv_3.txt")

rho1 = data[:, 0]
p1 = data[:, 1]
che1 = data[:, 2]
rho2 = data1[:, 0]
p2 = data1[:, 1]
che2 = data1[:, 2]
rho3 = data2[:, 0]
p3 = data2[:, 1]
che3 = data2[:, 2]

plt.figure(figsize=(10, 6))

plt.plot(rho1,p1, label="0")
plt.plot(rho2,p2, label ="100")
plt.plot(rho3,p3, label = "-100")
#plt.plot(rho1, che1, label="0")
#plt.plot(rho2, che2, label="100")
#plt.plot(rho3, che3, label="-100")
# plt.plot(rho_total,p5,label = "T")
plt.title("QCD Phase Diagram")
plt.xlabel("rho")

plt.legend()

plt.ylabel("Temperature (T)")
plt.grid(True)
plt.savefig("Pressure_rho.png")

plt.show()
