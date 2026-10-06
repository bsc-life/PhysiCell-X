"""Replay the 120-minute field through freshly built production transport code."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'transport_replay_120'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / 'tools'))
import compare_biological_outputs as comparison

stage = sys.argv[1]
assert stage in ['prepare', 'run']
harness = HERE.parent / 'reassessment_2026-09-09/transport_probe.cpp'
env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE', OPENBLAS_NUM_THREADS='1')

if stage == 'prepare':
    snapshot = comparison.read(HERE / 'runs/reference_coupled/output00000004.xml')
    positions = snapshot['labels']['position'].T
    bounds = snapshot['bounds']
    valid = np.all((positions >= bounds[:3]) & (positions < bounds[3:]), axis=1)
    cells = np.column_stack([positions, comparison.scalar(snapshot, 'total_volume'),
                             snapshot['labels']['net_export_rates'].T])[valid]
    np.savetxt(OUT / 'cells.txt', cells, fmt='%.17g')
    snapshot['micro'][[0, 1, 2, 4, 5, 6, 7, 8]].T.copy().tofile(OUT / 'initial.bin')
    metadata = dict(snapshot=str(snapshot['path']), time_min=120,
                    included_cells=int(valid.sum()), excluded_outside_cells=int((~valid).sum()),
                    steps=100, dt_min=.001, threads_per_process=1, builds=[])
    for name in ['x', 'reference']:
        build = HERE / ('build_' + name)
        records = json.loads((build / 'compile_initial.json').read_text())
        command = next(r['command'][:] for r in records if r['command'][-1].endswith('/test-cell-sources'))
        command = [str(harness) if arg.endswith('/test_cell_sources.cpp') else arg
                   for arg in command if arg != '-DSOURCE_TEST_MPI']
        if name == 'x': command.insert(1, '-DPROBE_MPI')
        with tempfile.TemporaryDirectory(prefix='late-transport-link-') as temporary:
            binary = Path(temporary) / ('probe_' + name)
            command[-1] = str(binary)
            result = subprocess.run(command, env=env, capture_output=True, text=True)
            record = dict(project=name, command=command, returncode=result.returncode,
                          stdout=result.stdout, stderr=result.stderr)
            metadata['builds'].append(record)
            (OUT / 'preparation.json').write_text(json.dumps(metadata, indent=2) + '\n')
            assert result.returncode == 0, result.stderr
            shutil.copy2(binary, OUT / ('probe_' + name))
    files = [harness, OUT / 'initial.bin', OUT / 'cells.txt', OUT / 'probe_x', OUT / 'probe_reference',
             HERE / 'build_x/provenance.json', HERE / 'build_reference/provenance.json'] + snapshot['files']
    (OUT / 'provenance.json').write_text(json.dumps({str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                                  for p in files}, indent=2) + '\n')
    print('Prepared 120-minute transport replay:', metadata['included_cells'], 'fixed cells', flush=True)
else:
    calls = []
    for sources in [0, 1]:
        for name, mode, ranks in [('reference', 'serial', 1), ('x', 'serial', 1), ('x', 'blocking', 2)]:
            tag = f'{name}_{mode}_sources{sources}'
            command = [str(OUT / ('probe_' + name)), mode, str(OUT / 'initial.bin'),
                       str(OUT / 'cells.txt'), str(OUT / tag), '1', '100', str(sources), '-1']
            if ranks == 2:
                command = ['mpirun', '--allow-run-as-root', '--oversubscribe', '--bind-to', 'none', '-n', '2'] + command
            result = subprocess.run(command, env=env, text=True, capture_output=True, timeout=120)
            calls.append(dict(tag=tag, command=command, returncode=result.returncode,
                              stdout=result.stdout, stderr=result.stderr))
            (OUT / 'calls.json').write_text(json.dumps(calls, indent=2) + '\n')
            assert result.returncode == 0, result.stderr
            print('Finished', tag, flush=True)
    def read(prefix):
        data = np.concatenate([np.fromfile(path).reshape(-1, 8) for path in sorted(OUT.glob(prefix + '_*.bin'))])
        data = data[np.lexsort((data[:, 0], data[:, 1], data[:, 2]))]
        assert data.shape == (250000, 8) and np.isfinite(data).all()
        return data
    results = []
    for sources in [0, 1]:
        reference = read('reference_serial_sources' + str(sources))
        for mode in ['serial', 'blocking']:
            actual = read('x_' + mode + '_sources' + str(sources))
            assert np.array_equal(actual[:, :3], reference[:, :3])
            error = actual[:, 3:] - reference[:, 3:]
            results.append(dict(mode=mode, fixed_sources=bool(sources), steps=100, dt_min=.001,
                                bitwise_equal=bool(np.array_equal(actual, reference)),
                                max_absolute_error=float(np.max(np.abs(error))),
                                max_error_by_substrate=np.max(np.abs(error), axis=0).tolist()))
    (OUT / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))
