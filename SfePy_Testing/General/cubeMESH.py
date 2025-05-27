import gmsh
import math
import meshio

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("CubeGroundMesh")

# Cube dimensions
cube_length = 2
cube_half_length = cube_length / 2

# Ground dimensions (larger than the cube)
ground_length = 6
ground_half_length = ground_length / 2

# Add the ground (large box) and the cube
ground = gmsh.model.occ.addBox(-ground_half_length, -ground_half_length, -ground_half_length, 
                               ground_length, ground_length, ground_length)
cube = gmsh.model.occ.addBox(-cube_half_length, -cube_half_length, -cube_half_length, 
                             cube_length, cube_length, cube_length)

# Subtract the cube from the ground to define their interface
gmsh.model.occ.cut([(3, ground)], [(3, cube)])
gmsh.model.occ.synchronize()

# Assign physical groups for the cube and ground
gmsh.model.addPhysicalGroup(3, [cube], name="Cube")   # The cold material (ice/water)
gmsh.model.addPhysicalGroup(3, [ground], name="Ground")  # The surrounding soil

# Get surface entities for boundary tagging
surfaces = gmsh.model.getEntities(dim=2)

# Tag boundaries for cube
for surf in surfaces:
    _, tag = surf
    com = gmsh.model.occ.getCenterOfMass(2, tag)
    z = com[2]
    if math.isclose(z, cube_half_length, rel_tol=1e-6):
        gmsh.model.addPhysicalGroup(2, [tag], name="Top")  # Top face of the cube
    elif math.isclose(z, -cube_half_length, rel_tol=1e-6):
        gmsh.model.addPhysicalGroup(2, [tag], name="Bottom")  # Bottom face of the cube
    elif abs(z) < cube_half_length:
        gmsh.model.addPhysicalGroup(2, [tag], name="Side")  # Side faces of the cube

# Generate the 3D mesh
gmsh.model.mesh.generate(3)

# Write the mesh to files
gmsh.write("cube_ground.msh")  # Gmsh MSH format
gmsh.write("cube_ground.vtk")  # VTK format

# Finalize Gmsh
gmsh.finalize()

# Convert the Gmsh MSH file to Dolfin-compatible formats using Meshio
mesh = meshio.read("cube_ground.msh")

# Save to Dolfin XML format
meshio.write("cube_ground.xml", mesh, file_format="dolfin-xml")

# Save to XDMF format (preferred for modern FEniCS/SfePy workflows)
meshio.write("cube_ground.xdmf", mesh, file_format="xdmf")
