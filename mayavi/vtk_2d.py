from mayavi import mlab
from mayavi.modules.surface import Surface

file_name = "solution-012.vtk"  # minimal example vtk file

# create a new figure, grab the engine that's created with it
#fig = mlab.figure()
engine = mlab.get_engine()
scene = engine.new_scene()

# open the vtk file, let mayavi figure it all out
vtk_file_reader = engine.open(file_name)

# plot surface corresponding to the data
surface = Surface()
engine.add_filter(surface, vtk_file_reader)

mlab.scalarbar(orientation='vertical')
mlab.axes(y_axis_visibility=False)

# Move the camera
scene.scene.camera.elevation(0)

# Save scene to image file
#mlab.savefig('image.pdf')

# block until figure is closed
mlab.show()
