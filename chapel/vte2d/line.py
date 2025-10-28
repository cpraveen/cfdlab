import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt

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
v = alongx.point_data['v']
x = alongx.point_data['Distance']

u = alongy.point_data['u']
y = alongy.point_data['Distance']

plt.figure(figsize=(10,5))

plt.subplot(121)
plt.plot(x, v, label="VTE: v(x,0.5)")
plt.plot(dv[:,0], dv[:,1], 'o', label="Ghia")
plt.xlabel("x"); plt.ylabel("v")
plt.legend()

plt.subplot(122)
plt.plot(u, y, label="VTE: u(0.5,y)")
plt.plot(du[:,1], du[:,0], 'o', label="Ghia")
plt.xlabel("u"); plt.ylabel("y")
plt.legend()

plt.savefig("vel.svg")
plt.show()
