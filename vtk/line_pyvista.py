import pyvista as pv

filename = "rect.vtk"

data = pv.read(filename)

# Define the line along which to sample the solution
# Example: line from point A to point B
point_a = [0.75, 0, 0]
point_b = [0.75, 1, 0]
n_points = 100  # Number of sample points along the line

data.plot_over_line(point_a, point_b, resolution=n_points, 
                    scalars="density",
                    title="Along line x=0.75",
                    ylabel="Density",
                    figsize=(10,5))

# Sample the mesh along the line
# TODO: How to get the data from this ?
sampled = data.sample_over_line(point_a, point_b, n_points)
print(sampled)
