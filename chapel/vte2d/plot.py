import pyvista as pv
from pylab import *

filename = "sol.vtk"

reader = pv.get_reader(filename)
data = reader.read()

#data.plot(scalars="psi", cpos='xy', window_size=(1200,1200), show_edges=False)

scalar = "psi"
cmap = "coolwarm"

psi = data.point_data["psi"]
pmin = psi.min();
pmax = psi.max();
print("psi min/max = ",pmin,pmax)

p = pv.Plotter(window_size=(2000, 2000))

c1 = linspace(pmin, 0, 20)
contours = data.contour(isosurfaces=c1, scalars=scalar)
p.add_mesh(contours, color="black", line_width=1.5,cmap=cmap)

c2 = linspace(0, pmax, 5)
contours = data.contour(isosurfaces=c2, scalars=scalar)
p.add_mesh(contours, color="black", line_width=1.5,cmap=cmap)

p.add_title("Stream function")
p.show_bounds(all_edges=True,
              location="outer",
              xtitle="x",
              ytitle="y",
              font_size=24,
              color=(0, 0, 0, 0.5))

p.view_xy()
p.show()
