import pyvista as pv

filename = "struct.vtk"

reader = pv.get_reader(filename)
data = reader.read()

# Alternately, you can read like this
# data = pv.read(filename)

data.plot(cpos='xy', window_size=(1200,1200), show_edges=True)
