import pyvista as pv
from pylab import *

filename = "sol.vtk"

with open(filename,"r") as f:
    f.readline()
    problem = f.readline().strip()

reader = pv.get_reader(filename)
data = reader.read()
print(data)

scalar = "psi"
cmap = "coolwarm"

psi = data.point_data["psi"]
pmin = psi.min();
pmax = psi.max();
print("psi min/max = ",pmin,pmax)
assert pmin < 0 and 0 < pmax, "Need pmin < 0 < pmax"

p = pv.Plotter(window_size=(1024, 1024))

# Contour levels of streamfunction from Ghia paper
cvals = [-1e-10, -1e-7, -1e-5, -1e-4, -0.01, -0.03, -0.05, -0.07, -0.09, \
         -0.1, -0.11, -0.115, -0.1175, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 2.5e-4, \
         5e-4, 1e-3, 1.5e-3, 3e-3]

contours = data.contour(isosurfaces=cvals, scalars=scalar)
p.add_mesh(contours, color="black", line_width=1.5,cmap=cmap)

p.add_title(problem + ", Streamlines")
p.show_bounds(all_edges=True,
              location="outer",
              xtitle="x",
              ytitle="y",
              font_size=28,
              color=(0, 0, 0, 0.5))

p.view_xy()
p.camera.zoom(1.1)
p.save_graphic("stream.svg") # pdf, eps, ps, svg
p.screenshot("stream.png")   # png, jpeg, jpg, bmp, tif, tiff
p.show()
