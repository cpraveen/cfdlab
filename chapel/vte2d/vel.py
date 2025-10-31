import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-Re", "--Re", type=int, help="Reynolds no.", required=True)
args = parser.parse_args()

Re = { 100: 1, 400: 2, 1000: 3, 3200: 4, 5000: 5, 7500: 6, 10000: 7 }
col = Re[args.Re]

filename = "sol.vtk"
data = pv.read(filename)

du = np.loadtxt("ghia_y_u.txt")
dv = np.loadtxt("ghia_x_v.txt")

# Define the line along which to sample the solution
# Example: line from point A to point B
p1 = [0.0, 0.5, 0.0]
p2 = [1.0, 0.5, 0.0]

p3 = [0.5, 0, 0]
p4 = [0.5, 1, 0]

npoints = 200  # Number of sample points along the line

# Sample the mesh along the line
alongx = data.sample_over_line(p1, p2, resolution=npoints)
alongy = data.sample_over_line(p3, p4, resolution=npoints)

# Plot using matplotlib
v = alongx.point_data["v"]
x = alongx.point_data["Distance"]

u = alongy.point_data["u"]
y = alongy.point_data["Distance"]

fig = plt.figure(figsize=(10,5))

plt.subplot(121)
plt.plot(x, v, label="VTE: v(x,0.5)")
plt.plot(dv[:,0], dv[:,col], "o", label="Ghia et al.")
plt.xlabel("x"); plt.ylabel("v")
plt.grid(True)
plt.legend()

plt.subplot(122)
plt.plot(u, y, label="VTE: u(0.5,y)")
plt.plot(du[:,col], du[:,0], "o", label="Ghia et al.")
plt.xlabel("u"); plt.ylabel("y")
plt.grid(True)
plt.legend()

fig.suptitle("Lid driven cavity, Re = " + str(args.Re))

plt.savefig("vel.svg")
plt.show()
