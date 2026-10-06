"""Rebuild/replay the focused tests. MPI requires permission to open local sockets.

Run from any directory with Python 3. Requires NumPy, SciPy, Matplotlib, g++, mpic++,
mpirun, and the two sibling projects. Does not modify their simulation outputs.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'transport_diagnosis'
sys.path.insert(0, str(ROOT / 'tools'))
import compare_biological_outputs as comparison


def main():
    calls = []
    env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE')
    def run(cmd):
        p = subprocess.run([str(v) for v in cmd], cwd=ROOT, env=env,
                           text=True, capture_output=True)
        calls.append(dict(command=[str(v) for v in cmd], returncode=p.returncode,
                          stdout=p.stdout, stderr=p.stderr))
        (OUT / 'reproduction_calls.json').write_text(json.dumps(calls, indent=2) + '\n')
        print(p.stdout, flush=True)
        if p.returncode:
            raise RuntimeError(p.stderr)

    OUT.mkdir(exist_ok=True)
    snap = comparison.read(ROOT.parent / 'PhysiCelldFBA/output/output00000001.xml')
    cells = np.column_stack([snap['labels']['position'].T, comparison.scalar(snap, 'total_volume'),
                             snap['labels']['net_export_rates'].T])
    np.savetxt(OUT / 'cells30.txt', cells, fmt='%.17g')
    initial = '/tmp/biological-transport-initial.bin'
    snap['micro'][[0, 1, 2, 4, 5, 6, 7, 8]].T.copy().tofile(initial)
    for project in ['reference', 'x']:
        run([sys.executable, HERE / 'build_transport_probe.py', project])
        binary = '/tmp/biological-transport-probe-' + project
        for threads, locked, partition in [(1, 0, -1), (2, 0, -1), (8, 0, -1),
                                           (112, 0, -1), (112, 1, -1)] + (
                                           [(2, 0, 0), (2, 0, 1)] if project == 'x' else []):
            tag = f'source_{project}_t{threads}_lock{locked}_part{partition}'
            run([binary, 'source', initial, OUT / 'cells30.txt', OUT / tag,
                 threads, 100, locked, partition])
        for replay in [0, 1]:
            run([binary, 'serial', initial, OUT / 'cells30.txt',
                 f'/tmp/biological-transport-{project}-replay{replay}', 1, 100, replay, -1])
    for replay in [0, 1]:
        prefix = '/tmp/biological-transport-blocking2' + ('-replay1' if replay else '')
        run(['mpirun', '--allow-run-as-root', '--oversubscribe', '-n', 2,
             '/tmp/biological-transport-probe-x', 'blocking', initial,
             OUT / 'cells30.txt', prefix, 1, 100, replay, -1])
    # Cache probe recompiles Basic_Agent, linking the existing supporting BioFVM objects.
    for project in ['reference', 'x']:
        location = ROOT if project == 'x' else ROOT.parent / 'PhysiCelldFBA'
        objects = ['BioFVM_vector.o', 'BioFVM_mesh.o', 'BioFVM_microenvironment.o',
                   'BioFVM_solvers.o', 'BioFVM_agent_container.o', 'BioFVM_utilities.o',
                   'BioFVM_matlab.o', 'pugixml.o']
        cmd = ['mpic++' if project == 'x' else 'g++', '-no-pie', '-std=c++17', '-O1',
               '-fopenmp', '-ffunction-sections', '-fdata-sections', '-Wl,--gc-sections',
               '-I' + str(location), HERE / 'cache_timestep_probe.cpp',
               location / 'BioFVM/BioFVM_basic_agent.cpp']
        if project == 'x':
            cmd += ['-DPROBE_MPI']
            objects += ['DistPhy_Cartesian.o', 'DistPhy_Environment.o']
        binary = '/tmp/biological-cache-probe-' + project
        run(cmd + [location / p for p in objects] + ['-o', binary])
        run([binary])
    run([sys.executable, HERE / 'analyse_transport_probes.py'])


if __name__ == '__main__':
    main()
