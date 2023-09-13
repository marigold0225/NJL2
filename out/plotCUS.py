import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("spinodal_G_sv_1.txt", delimiter=" ")
data1 = np.loadtxt("spinodal_G_sv_2.txt", delimiter=" ")
data2 = np.loadtxt("spinodal_G_sv_3.txt", delimiter=" ")

T1 = data[:, 1]
rho1 = data[:, 0]

T2 = data1[:, 1]
rho2 = data1[:, 0]

T3 = data2[:, 1]
rho3 = data2[:, 0]


plt.figure(figsize=(10, 6))

# plt.plot(rho1,p1, label="0")
# plt.plot(rho2,p2, label ="100")
# plt.plot(rho3,p3, label = "-100")
plt.scatter(x=rho1, y=T1, label="0")
plt.scatter(x=rho2, y=T2, label="100")
plt.scatter(x=rho3, y=T3, label="-200")

# plt.plot(rho_total,p5,label = "T")
plt.xlim(0, 1.5)
plt.title("QCD Phase Diagram")
plt.xlabel("rho")

plt.legend()

plt.ylabel("Temperature (T)")
plt.grid(True)
plt.savefig("spin_C++.png")

plt.show()
