#!/usr/bin/env python3
"""Compile and run an isolated mesh regression; no simulation or MPI launch."""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]

with tempfile.TemporaryDirectory(prefix="physicell-voxel-index-") as directory:
    binary = Path(directory) / "voxel_local_index"
    subprocess.run([
        "mpic++", "-std=c++11", "-O1", "-fopenmp",
        "-ffunction-sections", "-fdata-sections", "-Wl,--gc-sections",
        "tests/unit/voxel_local_index.cpp", "BioFVM/BioFVM_mesh.cpp",
        "BioFVM/BioFVM_vector.cpp", "-o", str(binary),
    ], cwd=ROOT, check=True)
    subprocess.run([str(binary)], cwd=ROOT, check=True)
