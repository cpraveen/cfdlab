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

# negative contours
c1 = linspace(pmin, 0, 20)
contours = data.contour(isosurfaces=c1, scalars=scalar)
p.add_mesh(contours, color="black", line_width=1.5,cmap=cmap)

# positive contours, these show corner vortices
c2 = linspace(0, pmax, 5)
contours = data.contour(isosurfaces=c2, scalars=scalar)
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
p.save_graphic("psi.svg") # pdf, eps, ps, svg
p.screenshot("psi.png")   # png, jpeg, jpg, bmp, tif, tiff
p.show()
