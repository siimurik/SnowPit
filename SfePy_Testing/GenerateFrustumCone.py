import gmsh
import meshio
import math as m

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("FrustumCone")

# Parameters for the frustum
bottom_radius = 1.0
top_radius = 4.0
height = 3.0
num_segments = 36  # Number of segments for the circle

# Helper function to add circle points
def add_circle_points(radius, z, num_segments):
    points = []
    for i in range(num_segments):
        angle = 2 * m.pi * i / num_segments
        x = radius * m.cos(angle)
        y = radius * m.sin(angle)
        points.append(gmsh.model.occ.addPoint(x, y, z))
    return points

# Add points for the bottom and top circles
bottom_points = add_circle_points(bottom_radius, 0.0, num_segments)
top_points = add_circle_points(top_radius, height, num_segments)

# Add lines for the bottom and top circles
bottom_lines = [
    gmsh.model.occ.addLine(bottom_points[i], bottom_points[(i + 1) % num_segments])
    for i in range(num_segments)
]
top_lines = [
    gmsh.model.occ.addLine(top_points[i], top_points[(i + 1) % num_segments])
    for i in range(num_segments)
]

# Add vertical lines connecting the top and bottom circles
vertical_lines = [
    gmsh.model.occ.addLine(bottom_points[i], top_points[i]) for i in range(num_segments)
]

# Create curve loops for the top and bottom
bottom_loop = gmsh.model.occ.addCurveLoop(bottom_lines)
top_loop = gmsh.model.occ.addCurveLoop(top_lines)

# Add surfaces for the top and bottom
bottom_surface = gmsh.model.occ.addPlaneSurface([bottom_loop])
top_surface = gmsh.model.occ.addPlaneSurface([top_loop])

# Add side surfaces
side_surfaces = []
for i in range(num_segments):
    side_loop = gmsh.model.occ.addCurveLoop([
        bottom_lines[i],
        vertical_lines[(i + 1) % num_segments],
        -top_lines[i],
        -vertical_lines[i],
    ])
    side_surfaces.append(gmsh.model.occ.addPlaneSurface([side_loop]))

# Create a surface loop and volume
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
output_file = "frustum_cone.mesh"
meshio.write(output_file, mesh, file_format="medit")

# Finalize Gmsh
gmsh.finalize()

print(f"Tetrahedral mesh written to '{output_file}'")
