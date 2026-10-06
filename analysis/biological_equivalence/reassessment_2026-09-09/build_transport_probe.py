"""Build a focused probe from either project's production BioFVM source."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
name = sys.argv[1]
assert name in ['x', 'reference']
project = ROOT if name == 'x' else ROOT.parent / 'PhysiCelldFBA'
sources = [project / 'BioFVM' / s for s in [
    'BioFVM_vector.cpp', 'BioFVM_mesh.cpp', 'BioFVM_microenvironment.cpp',
    'BioFVM_solvers.cpp', 'BioFVM_basic_agent.cpp', 'BioFVM_agent_container.cpp',
    'BioFVM_utilities.cpp', 'BioFVM_matlab.cpp']]
cmd = ['mpic++' if name == 'x' else 'g++', '-std=c++17', '-O2', '-fopenmp',
       '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections', '-mavx2', '-mfma',
       '-I' + str(project), str(HERE / 'transport_probe.cpp')]
if name == 'x':
    cmd += ['-DPROBE_MPI']
    sources += [project / 'DistPhy' / s for s in ['DistPhy_Cartesian.cpp', 'DistPhy_Environment.cpp']]
cmd += [str(p) for p in sources] + ['-o', '/tmp/biological-transport-probe-' + name]
result = subprocess.run(cmd, text=True, capture_output=True)
record = dict(command=cmd, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
              source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sources + [HERE / 'transport_probe.cpp']})
(HERE / 'transport_diagnosis' / ('build_' + name + '.json')).write_text(json.dumps(record, indent=2) + '\n')
print(name, result.returncode, result.stdout, result.stderr[-5000:])
sys.exit(result.returncode)
