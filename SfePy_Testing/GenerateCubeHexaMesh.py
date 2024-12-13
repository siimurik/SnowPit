import numpy as np

# Parameters for the cube
cube_size = 1
n_divisions = 10  # Number of divisions per edge

# Derived parameters
step = cube_size / n_divisions
vertices = []
hexahedra = []

# Generate vertices
vertex_index = 1
vertex_map = {}
for i in range(n_divisions + 1):
    for j in range(n_divisions + 1):
        for k in range(n_divisions + 1):
            x = -0.5 * cube_size + i * step
            y = -0.5 * cube_size + j * step
            z = -0.5 * cube_size + k * step
            vertices.append(f"{x:.10e} {y:.10e} {z:.10e} 0")
            vertex_map[(i, j, k)] = vertex_index
            vertex_index += 1

# Generate hexahedra
for i in range(n_divisions):
    for j in range(n_divisions):
        for k in range(n_divisions):
            # Define the 8 corners of the hexahedron
            v0 = vertex_map[(i, j, k)]
            v1 = vertex_map[(i + 1, j, k)]
            v2 = vertex_map[(i + 1, j + 1, k)]
            v3 = vertex_map[(i, j + 1, k)]
            v4 = vertex_map[(i, j, k + 1)]
            v5 = vertex_map[(i + 1, j, k + 1)]
            v6 = vertex_map[(i + 1, j + 1, k + 1)]
            v7 = vertex_map[(i, j + 1, k + 1)]

            hexahedra.append(f"{v0} {v1} {v2} {v3} {v4} {v5} {v6} {v7} 0")

# Write to .mesh file
with open("cubeHexa.mesh", "w") as f:
    f.write("MeshVersionFormatted 2\n")
    f.write("Dimension 3\n")

    # Write vertices
    f.write("Vertices\n")
    f.write(f"{len(vertices)}\n")
    f.write("\n".join(vertices) + "\n")

    # Write hexahedra
    f.write("Hexahedra\n")
    f.write(f"{len(hexahedra)}\n")
    f.write("\n".join(hexahedra) + "\n")

print("Mesh file 'cubeHexa.mesh' created successfully.")