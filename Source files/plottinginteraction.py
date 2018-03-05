import numpy as np
import matplotlib.pyplot as plt

eigenvec1 = []
eigenvec2 = []
with open("eigenvecs_noninteract.txt","r") as data:
    for line in data:
        numbers_float = list(map(float, line.split()))
        # print(numbers_float)
        eigenvec1.append(numbers_float[0])
with open("eigenvecs_interact.txt","r") as data:
    for line in data:
        numbers_float = list(map(float, line.split()))
        # print(numbers_float)
        eigenvec2.append(numbers_float[0])
rho = []
with open("rhovals.txt","r") as data:
    for line in data:
        #numbers_float = list(map(float, line.split()))
        rho.append(float(line))



eigenvec1 = np.asarray(eigenvec1)
eigenvec2 = np.asarray(eigenvec2)
rho = np.asarray(rho)

eigenvec1 = (eigenvec1*rho)**2
eigenvec2 = (eigenvec2*rho)**2

area1 = np.trapz(eigenvec1, dx=(rho[1]-rho[0]))
area2 = np.trapz(eigenvec2, dx=(rho[1]-rho[0]))

eigenvec1 = eigenvec1 / area1
eigenvec2 = eigenvec2 / area2

area1 = np.trapz(eigenvec1, dx=(rho[1]-rho[0]))
area2 = np.trapz(eigenvec2, dx=(rho[1]-rho[0]))
print("Area of wavefuntion: ", area1, area2)

plt.plot(rho,eigenvec1, label="non-interacting")
plt.plot(rho,eigenvec2, label="interacting")
plt.legend()
plt.ylabel("wavfunction r^2*R^2")
plt.xlabel("radial distance r")
plt.title("probabilities for the interacting and non-interacting, w = 5")
#
plt.savefig("intervsnon_w5.png")
plt.show()
