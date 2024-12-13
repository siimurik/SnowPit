import meshio

# Input and output file paths
input_file = "cubeTetra.mesh"
output_file = "cubeTetra_converted.mesh"

# Read Gmsh format
mesh = meshio.read(input_file)

# Write Medit-compatible .mesh format
meshio.write(output_file, mesh, file_format="medit")
print(f"Converted mesh written to '{output_file}'")
