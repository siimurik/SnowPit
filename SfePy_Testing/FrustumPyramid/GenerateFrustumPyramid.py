import gmsh
import meshio

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("TruncatedPyramid")

# Parameters for the shape
top_length = 10.0
top_width  = 16.0
bottom_length = 4.0
bottom_width  = 8.0
height = 4.0
char_length = 0.90  # Characteristic length for mesh refinement

# Add points for the top rectangle
top_points = [
    gmsh.model.occ.addPoint(-top_length / 2, -top_width / 2, height, char_length),
    gmsh.model.occ.addPoint(top_length / 2, -top_width / 2, height, char_length),
    gmsh.model.occ.addPoint(top_length / 2, top_width / 2, height, char_length),
    gmsh.model.occ.addPoint(-top_length / 2, top_width / 2, height, char_length),
]

# Add points for the bottom rectangle
bottom_points = [
    gmsh.model.occ.addPoint(-bottom_length / 2, -bottom_width / 2, 0, char_length),
    gmsh.model.occ.addPoint(bottom_length / 2, -bottom_width / 2, 0, char_length),
    gmsh.model.occ.addPoint(bottom_length / 2, bottom_width / 2, 0, char_length),
    gmsh.model.occ.addPoint(-bottom_length / 2, bottom_width / 2, 0, char_length),
]

# Add lines for the top and bottom rectangles
top_lines = [
    gmsh.model.occ.addLine(top_points[i], top_points[(i + 1) % 4]) for i in range(4)
]
bottom_lines = [
    gmsh.model.occ.addLine(bottom_points[i], bottom_points[(i + 1) % 4]) for i in range(4)
]

# Add vertical lines connecting top and bottom rectangles
vertical_lines = [
    gmsh.model.occ.addLine(bottom_points[i], top_points[i]) for i in range(4)
]

# Add curve loops and surfaces for the top and bottom rectangles
top_loop = gmsh.model.occ.addCurveLoop(top_lines)
bottom_loop = gmsh.model.occ.addCurveLoop(bottom_lines)

# Add surfaces for the top and bottom
top_surface = gmsh.model.occ.addPlaneSurface([top_loop])
bottom_surface = gmsh.model.occ.addPlaneSurface([bottom_loop])

# Create a surface loop and volume
side_surfaces = [
    gmsh.model.occ.addPlaneSurface(
        [gmsh.model.occ.addCurveLoop([bottom_lines[i], vertical_lines[i], -top_lines[i], -vertical_lines[(i + 1) % 4]])]
    )
    for i in range(4)
]
surface_loop = gmsh.model.occ.addSurfaceLoop([top_surface, bottom_surface] + side_surfaces)
volume = gmsh.model.occ.addVolume([surface_loop])

# Synchronize and generate the mesh
gmsh.model.occ.synchronize()
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
output_file = "frustum_pyramid.mesh"
meshio.write(output_file, mesh, file_format="medit")

# Finalize Gmsh
gmsh.finalize()

print(f"Tetrahedral mesh written to '{output_file}'")