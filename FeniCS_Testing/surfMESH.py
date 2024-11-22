import gmsh
import math
import meshio

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("CubeMesh")

# Define a cube geometry (centered at the origin with side length 2)
length = 2
half_length = length / 2
gmsh.model.occ.addBox(-half_length, -half_length, -half_length, length, length, length)
gmsh.model.occ.synchronize()

# Define physical groups for boundaries (useful for boundary conditions)
surfaces = gmsh.model.getEntities(dim=2)

for surf in surfaces:
    _, tag = surf
    com = gmsh.model.occ.getCenterOfMass(2, tag)
    z = com[2]
    if math.isclose(z, half_length, rel_tol=1e-6):
        gmsh.model.addPhysicalGroup(2, [tag], name="Top")
    elif math.isclose(z, -half_length, rel_tol=1e-6):
        gmsh.model.addPhysicalGroup(2, [tag], name="Bottom")
    else:
        gmsh.model.addPhysicalGroup(2, [tag], name="Side")

# Generate the 3D mesh (tetrahedral)
gmsh.model.mesh.generate(3)

# Write the mesh to a Gmsh MSH file
gmsh.write("outputTetMesh.msh")  # Gmsh MSH format
gmsh.write("outputTetMesh.vtk")  # VTK format

# Finalize Gmsh
gmsh.finalize()

# Convert the .msh file to Dolfin XML format using meshio
mesh = meshio.read("outputTetMesh.msh")
#meshio.write("outputTetMesh.xml", mesh, file_format="dolfin-xml")
#cmeshio.write("outputTetMesh.xdmf", mesh, file_format="xdmf")

