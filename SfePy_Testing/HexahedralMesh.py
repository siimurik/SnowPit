import gmsh
import meshio

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("CubeHexaMesh")

# Create a cube geometry with dimensions 1x1x1 centered at the origin
cube_size = 1.0
half_size = cube_size / 2.0
gmsh.model.occ.addBox(-half_size, -half_size, -half_size, cube_size, cube_size, cube_size)
gmsh.model.occ.synchronize()

# Apply transfinite meshing to the cube to create structured hexahedral elements
box_volume = gmsh.model.getEntities(dim=3)[0][1]  # Get the cube's volume tag
gmsh.model.mesh.setTransfiniteVolume(box_volume)
gmsh.model.mesh.setRecombine(3, box_volume)  # Enable hexahedral elements

# Generate the hexahedral mesh
gmsh.model.mesh.generate(3)

# Retrieve mesh data directly from Gmsh
nodes = gmsh.model.mesh.getNodes()
elements = gmsh.model.mesh.getElements(dim=3)

# Extract node coordinates and hexahedral elements
node_coords = nodes[1].reshape((-1, 3))
hex_elements = elements[2][0] - 1  # Convert 1-based indexing to 0-based for Python

# Create a Meshio mesh object
mesh = meshio.Mesh(
    points=node_coords,
    cells=[("hexahedron", hex_elements.reshape((-1, 8)))]  # 8 nodes per hexahedron
)

# Write the final mesh
