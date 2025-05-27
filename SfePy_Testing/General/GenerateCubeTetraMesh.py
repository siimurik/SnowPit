import gmsh
import meshio

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("CubeTetraMesh")

# Create a cube geometry with dimensions 1x1x1 centered at the origin
cube_size = 1.0
half_size = cube_size / 2.0
gmsh.model.occ.addBox(-half_size, -half_size, -half_size, cube_size, cube_size, cube_size)
gmsh.model.occ.synchronize()

# Generate the tetrahedral mesh
gmsh.model.mesh.generate(3)

# Retrieve mesh data directly from Gmsh
nodes = gmsh.model.mesh.getNodes()
elements = gmsh.model.mesh.getElements(dim=3)

# Extract node coordinates and tetrahedral elements
node_coords = nodes[1].reshape((-1, 3))
tetra_elements = elements[2][0] - 1  # Convert 1-based indexing to 0-based for Python

# Create a Meshio mesh object
mesh = meshio.Mesh(
    points=node_coords,
    cells=[("tetra", tetra_elements.reshape((-1, 4)))]
)

# Write the final mesh directly to Medit format
output_file = "cubeTetra_converted.mesh"
meshio.write(output_file, mesh, file_format="medit")

# Finalize Gmsh
gmsh.finalize()

print(f"Converted tetrahedral mesh written to '{output_file}'")
