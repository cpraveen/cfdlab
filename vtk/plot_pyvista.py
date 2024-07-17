import pyvista as pv

filename = "struct.vtk"
reader = pv.get_reader(filename)
data = reader.read()
data.plot(cpos='xy', window_size=(1200,1200), show_edges=True)
