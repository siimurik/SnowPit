from mpi4py import MPI
from dolfinx import mesh
import numpy as np

# Create a simple unit square mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)
print(f"Mesh created with {domain.topology.index_map(0).size_local} cells.")
