"""Test whether a paired spatial residual follows the saved source layout.

Usage: replay_source_layout.py [snapshot minute, default 120] [--reuse-results]
The reference transport binary is reused from the verified 120-minute control.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
time_min = int(next((a for a in sys.argv[1:] if a != '--reuse-results'), '120'))
assert time_min in [30, 60, 90, 120, 150, 180]
OUT = HERE / ('transport_replay_' + str(time_min))
OUT.mkdir(exist_ok=True)
binary = HERE / 'transport_replay_120/probe_reference'
sys.path.insert(0, str(ROOT / 'tools'))
import compare_biological_outputs as comparison

snapshots = {name: comparison.read(HERE / 'runs' / (name + '_coupled') / f'output{time_min // 30:08d}.xml')
             for name in ['x', 'reference']}
cells = {}
for name, snap in snapshots.items():
    positions = snap['labels']['position'].T
    valid = np.all((positions >= snap['bounds'][:3]) & (positions < snap['bounds'][3:]), axis=1)
    cells[name] = np.column_stack([positions, comparison.scalar(snap, 'total_volume'),
                                   snap['labels']['net_export_rates'].T])[valid]
specific_rate = cells['reference'][:, 4:].sum(axis=0) / cells['reference'][:, 3].sum()
metadata = dict(time_min=time_min, steps=100, dt_min=.001, solver='reference serial',
                common_rate_per_volume=specific_rate.tolist(),
                cells={name: len(data) for name, data in cells.items()},
                outside_cells={name: len(snapshots[name]['labels']['ID'][0]) - len(data) for name, data in cells.items()},
                binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(), calls=[])
env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE', OPENBLAS_NUM_THREADS='1')
if '--reuse-results' not in sys.argv[1:]:
    snapshots['reference']['micro'][[0, 1, 2, 4, 5, 6, 7, 8]].T.copy().tofile(OUT / 'initial.bin')
    metadata['input_sha256'] = hashlib.sha256((OUT / 'initial.bin').read_bytes()).hexdigest()
cases = [] if '--reuse-results' in sys.argv[1:] else [('reference', False), ('x', False), ('x', True), ('reference', True)]
for name, common_rate in cases:
    tag = 'layout_' + name + ('_common_rate' if common_rate else '_actual_rate')
    if name == 'reference' and not common_rate: tag = 'reference_serial_sources1'
    data = cells[name].copy()
    if common_rate:
        data[:, 4:] = data[:, 3, None] * specific_rate
    fixture = OUT / (tag + '.txt')
    np.savetxt(fixture, data, fmt='%.17g')
    command = [str(binary), 'serial', str(OUT / 'initial.bin'),
               str(fixture), str(OUT / tag), '1', '100', '1', '-1']
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=120)
    metadata['calls'].append(dict(command=command, returncode=result.returncode,
                                  stdout=result.stdout, stderr=result.stderr,
                                  fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest()))
    (OUT / 'layout_calls.json').write_text(json.dumps(metadata, indent=2) + '\n')
    assert result.returncode == 0, result.stderr


def read(name):
    data = np.fromfile(OUT / (name + '_0.bin')).reshape(-1, 8)
    data = data[np.lexsort((data[:, 0], data[:, 1], data[:, 2]))]
    assert data.shape == (250000, 8) and np.isfinite(data).all()
    assert np.array_equal(data[:, :3], snapshots['reference']['micro'][:3].T)
    return data[:, 3:]


observed = (snapshots['x']['micro'][4:] - snapshots['reference']['micro'][4:]).T
predictions = {
    'actual_saved_rates': read('layout_x_actual_rate') - read('reference_serial_sources1'),
    'common_rates_per_cell_volume': read('layout_x_common_rate') - read('layout_reference_common_rate')}
rows = []
for method, predicted in predictions.items():
    for d, substrate in enumerate(snapshots['reference']['names']):
        denominator = np.linalg.norm(observed[:, d])
        residual = predicted[:, d] - observed[:, d]
        rows.append(dict(method=method, substrate=substrate,
                         observed_difference_l2=float(denominator),
                         predicted_difference_l2=float(np.linalg.norm(predicted[:, d])),
                         residual_relative_l2=float(np.linalg.norm(residual) / denominator) if denominator else None,
                         max_residual=float(np.max(np.abs(residual)))))
(OUT / 'layout_summary.json').write_text(json.dumps(rows, indent=2) + '\n')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/physicell-comparison-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

axes_values = [np.unique(snapshots['reference']['micro'][d]) for d in range(3)]
shape = tuple(len(axis) for axis in axes_values[::-1]) + (5,)
observed_xy = observed.reshape(shape).mean(axis=0)
predicted_xy = predictions['common_rates_per_cell_volume'].reshape(shape).mean(axis=0)
window = np.abs(axes_values[0]) <= 25
fig, axes = plt.subplots(4, 3, figsize=(11, 11), constrained_layout=True)
for row, d in enumerate([0, 1, 2, 4]):
    actual = observed_xy[:, :, d][window][:, window]
    predicted = predicted_xy[:, :, d][window][:, window]
    bound = max(float(np.max(np.abs(actual))), 1e-16)
    for column, data in enumerate([actual, predicted, predicted - actual]):
        plot = axes[row, column].imshow(data, origin='lower', extent=(-26, 26, -26, 26),
                                        cmap='RdBu_r', vmin=-bound, vmax=bound)
        axes[row, column].set_xticks([-20, 0, 20])
        axes[row, column].set_yticks([-20, 0, 20])
        if row == 3: axes[row, column].set_xlabel('x (µm)')
    axes[row, 0].set_ylabel(snapshots['reference']['names'][d] + '\ny (µm)')
    fig.colorbar(plot, ax=axes[row, :], shrink=.8, label='Saved concentration units')
for ax, title in zip(axes[0], ['Observed coupled difference', 'Common-rate layout replay', 'Replay minus observed']):
    ax.set_title(title, fontsize=10)
fig.suptitle(f'{time_min}-minute concentration differences: PhysiCell-X − reference\n'
             'Mean over all four z layers; identical color scale within each row')
fig.savefig(OUT / 'layout_replay.png', dpi=180)
plt.close(fig)
print(json.dumps(rows, indent=2))
