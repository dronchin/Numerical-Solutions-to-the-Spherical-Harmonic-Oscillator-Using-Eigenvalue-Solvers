import numpy as np
import matplotlib.pyplot as plt

eigenvec1 = []
eigenvec2 = []
eigenvec3 = []
with open("eigenvecs_interact.txt","r") as data:
    for line in data:
        numbers_float = list(map(float, line.split()))
        # print(numbers_float)
        eigenvec1.append(numbers_float[0])
        eigenvec2.append(numbers_float[1])
        eigenvec3.append(numbers_float[2])
rho = []
with open("rhovals.txt","r") as data:
    for line in data:
        #numbers_float = list(map(float, line.split()))
        rho.append(float(line))


eigenvec1 = np.asarray(eigenvec1)
eigenvec2 = np.asarray(eigenvec2)
eigenvec3 = np.asarray(eigenvec3)
rho = np.asarray(rho)

eigenvec1 = (eigenvec1*rho)**2
eigenvec2 = (eigenvec2*rho)**2
eigenvec3 = (eigenvec3*rho)**2

area1 = np.trapz(eigenvec1, dx=(rho[1]-rho[0]))
area2 = np.trapz(eigenvec2, dx=(rho[1]-rho[0]))
area3 = np.trapz(eigenvec3, dx=(rho[1]-rho[0]))

eigenvec1 = eigenvec1 / area1
eigenvec2 = eigenvec2 / area2
eigenvec3 = eigenvec3 / area3

area1 = np.trapz(eigenvec1, dx=(rho[1]-rho[0]))
area2 = np.trapz(eigenvec2, dx=(rho[1]-rho[0]))
area3 = np.trapz(eigenvec3, dx=(rho[1]-rho[0]))

print("Area of each wavefunction: ", area1,area2,area3)

# print(len(rho),len(eigenvec1))
plt.plot(rho,eigenvec1)
plt.plot(rho,eigenvec2)
plt.plot(rho,eigenvec3)
plt.ylabel("wavfunction r^2*R^2")
plt.xlabel("radial distance r")
plt.title("Probabilities for the first three energies")

plt.savefig("3probabilities.png")
plt.show()
