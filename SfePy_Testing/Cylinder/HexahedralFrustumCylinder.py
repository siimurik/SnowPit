import gmsh
import meshio
import numpy as np

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("HexahedralFrustumCylinder")

# Parameters
bottom_radius = 1.0
top_radius = 2.0
height = 3.0
num_segments = 12
divisions = 4

# Create points
bottom_points, top_points = [], []
for i in range(num_segments):
    angle = 2 * np.pi * i / num_segments
    bottom_points.append(gmsh.model.occ.addPoint(bottom_radius * np.cos(angle), bottom_radius * np.sin(angle), 0))
    top_points.append(gmsh.model.occ.addPoint(top_radius * np.cos(angle), top_radius * np.sin(angle), height))

# Create lines
bottom_lines = [gmsh.model.occ.addLine(bottom_points[i], bottom_points[(i + 1) % num_segments]) for i in range(num_segments)]
top_lines = [gmsh.model.occ.addLine(top_points[i], top_points[(i + 1) % num_segments]) for i in range(num_segments)]

side_lines = []
side_surfaces = []
for i in range(num_segments):
    l1 = gmsh.model.occ.addLine(bottom_points[i], top_points[i])
    l2 = gmsh.model.occ.addLine(bottom_points[(i + 1) % num_segments], top_points[(i + 1) % num_segments])
    side_lines.append((l1, l2))
    side_loop = gmsh.model.occ.addCurveLoop([bottom_lines[i], l2, -top_lines[i], -l1])
    side_surfaces.append(gmsh.model.occ.addPlaneSurface([side_loop]))

# Create top and bottom surfaces
bottom_loop = gmsh.model.occ.addCurveLoop(bottom_lines)
top_loop = gmsh.model.occ.addCurveLoop(top_lines)
bottom_surface = gmsh.model.occ.addPlaneSurface([bottom_loop])
top_surface = gmsh.model.occ.addPlaneSurface([top_loop])

# Synchronize geometry
gmsh.model.occ.synchronize()

# Transfinite constraints
for line in bottom_lines + top_lines + [l for pair in side_lines for l in pair]:
    gmsh.model.mesh.setTransfiniteCurve(line, divisions)

gmsh.model.mesh.setTransfiniteSurface(bottom_surface, cornerTags=bottom_points)
gmsh.model.mesh.setTransfiniteSurface(top_surface, cornerTags=top_points)

for i, surface in enumerate(side_surfaces):
    gmsh.model.mesh.setTransfiniteSurface(
        surface,
        cornerTags=[
            bottom_points[i],
            bottom_points[(i + 1) % num_segments],
            top_points[(i + 1) % num_segments],
            top_points[i],
        ]
    )

# Create volume
surface_loop = gmsh.model.occ.addSurfaceLoop([bottom_surface, top_surface] + side_surfaces)
volume = gmsh.model.occ.addVolume([surface_loop])

gmsh.model.occ.synchronize()

# Apply transfinite volume and recombine for hexahedral mesh
gmsh.model.mesh.setTransfiniteVolume(volume)
gmsh.model.mesh.setRecombine(3, volume)

# Generate and save mesh
gmsh.model.mesh.generate(3)
gmsh.write("frustum_cylinder_hex.msh")

# Convert with meshio
mesh = meshio.read("frustum_cylinder_hex.msh")
meshio.write("frustum_cylinder_hex.mesh", mesh, file_format="medit")

gmsh.finalize()
print("Hexahedral mesh written to 'frustum_cylinder_hex.mesh'")
