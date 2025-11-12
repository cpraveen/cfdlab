import pyvista as pv
from pylab import *

filename = "sol.vtk"
reader = pv.get_reader(filename)
data = reader.read()
print(data)

scalar = "sol"
cmap = "coolwarm"

p = pv.Plotter(window_size=(1024, 1024))

contours = data.contour(isosurfaces=20, scalars=scalar)
p.add_mesh(contours, color="black", line_width=1.5, cmap=cmap)

p.show_bounds(all_edges=True,
              location="outer",
              xtitle="x",
              ytitle="y",
              font_size=24,
              color=(0, 0, 0, 0.5))

p.view_xy()
p.camera.zoom(1.1)
p.show()
