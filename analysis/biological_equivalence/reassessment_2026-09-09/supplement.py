"""Reproduce before/after tables, figures and focused validation from saved outputs."""
import csv
import json
import os
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/physicell-comparison-matplotlib')
HERE = Path(__file__).resolve().parent
OLD = HERE.parent
ROOT = OLD.parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
import compare_biological_outputs as c
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rows(folder, name):
    with (folder / name).open() as f:
        return list(csv.DictReader(f))


def index(data, *keys):
    return {tuple(r[k] for k in keys): r for r in data}


def settings(path):
    result = {}
    def visit(node, key):
        result[key] = [(node.text or '').strip(), node.attrib]
        for i, child in enumerate(node):
            visit(child, key + '/' + child.tag + '[' + str(i) + ']')
    visit(ET.parse(path).getroot(), 'root')
    return result


def main():
    old_cells = index(rows(OLD, 'cell_metrics.csv'), 'time_min', 'project')
    new_cells = index(rows(HERE, 'cell_metrics.csv'), 'time_min', 'project')
    old_fields = index(rows(OLD, 'substrate_metrics.csv'), 'time_min', 'substrate')
    new_fields = index(rows(HERE, 'substrate_metrics.csv'), 'time_min', 'substrate')
    before_after = []
    for (time, project), previous in old_cells.items():
        if project != 'PhysiCell-X':
            continue
        current = new_cells[time, project]
        reference = new_cells[time, 'PhysiCelldFBA']
        for metric in ['n_cells', 'total_volume', 'pressure_mean', 'pressure_max',
                       'radius_gyration', 'overlap_max', 'severe_overlap_pairs',
                       'severe_overlap_pairs_crossing_x_zero', 'geometric_pressure_mean']:
            before_after.append(dict(time_min=float(time), component='cells', metric=metric,
                                     before=float(previous[metric]), after=float(current[metric]),
                                     reference=float(reference[metric])))
    for (time, substrate), previous in old_fields.items():
        for metric in ['integrated_change_relative_difference', 'perturbation_relative_l2',
                       'gradient_relative_l2']:
            if previous[metric]:
                before_after.append(dict(time_min=float(time), component=substrate, metric=metric,
                                         before=float(previous[metric]),
                                         after=float(new_fields[time, substrate][metric]), reference=0.))
    c.csv_write(HERE / 'before_after.csv', before_after)
    provenance = json.loads((HERE / 'provenance.json').read_text())
    previous_provenance = json.loads((OLD / 'provenance.json').read_text())
    checks = {}
    for side in ['project', 'reference']:
        # The earlier analysis ran at a different absolute checkout path.
        files = [r for r in previous_provenance['files']
                 if str(Path(r['path']).parent) == previous_provenance[side]]
        changed = [Path(r['path']).name for r in files
                   if c.sha256(Path(provenance[side]) / Path(r['path']).name) != r['sha256']]
        checks[side + '_baseline_hash_comparison'] = dict(checked=len(files), changed=changed)
    assert not checks['reference_baseline_hash_comparison']['changed']
    xs = settings(Path(provenance['project']) / 'PhysiCell_settings.xml')
    rs = settings(Path(provenance['reference']) / 'PhysiCell_settings.xml')
    checks['settings_differences'] = {k: [xs.get(k), rs.get(k)] for k in sorted(xs.keys() | rs.keys())
                                      if xs.get(k) != rs.get(k)}
    for key, previous in old_cells.items():
        if key[1] == 'PhysiCelldFBA':
            assert previous == new_cells[key], key
    checks['reference_metrics_unchanged_at_original_times'] = True
    assert all(float(r['dead_count']) == 0 and float(r['velocity_max']) == 0
               for r in new_cells.values())
    assert all(json.loads(r['is_motile_counts']) == {'0': int(r['n_cells'])}
               and json.loads(r['current_phase_counts']) == {'14': int(r['n_cells'])}
               and json.loads(r['cycle_model_counts']) == {'5': int(r['n_cells'])}
               for r in new_cells.values())
    checks['all_cells_alive_nonmotile_phase14_cycle5_and_saved_velocities_zero'] = True
    assert all(float(r['x_negative_voxels']) == float(r['reference_negative_voxels']) == 0
               for r in new_fields.values())
    assert all(float(r['rmse']) == 0 for r in new_fields.values() if r['substrate'] == 'lactate')
    checks['nonnegative_fields_and_inactive_lactate'] = True
    final_time = max(provenance['common_times'])
    a = c.read(c.snapshots(Path(provenance['project']))[final_time])
    b = c.read(c.snapshots(Path(provenance['reference']))[final_time])
    baseline = c.read(c.snapshots(Path(provenance['reference']))[0])['micro'][4:]
    self_comparison = c.field_metrics(b, b, final_time, baseline)
    assert all(r['rmse'] == r['max_absolute_error'] == 0 and r['gradient_relative_l2'] in [0, None]
               for r in self_comparison)
    checks['same_input_final_field_concentration_and_gradient_errors_zero'] = True
    source_paths = ['BioFVM/BioFVM_mesh.cpp', 'core/PhysiCell_cell_container.cpp',
                    'core/PhysiCell_cell.cpp', 'core/PhysiCell_cell_container_mechanics.cpp',
                    'core/PhysiCell_standard_models.cpp', 'main.cpp', 'custom_modules/custom.cpp',
                    'config/Ecoli_core.xml', 'bacterial-colony', 'test-45570163.log',
                    'tools/compare_biological_outputs.py']
    checks['current_source_binary_log_hashes'] = {p: c.sha256(ROOT / p) for p in source_paths}
    checks['reference_model_hash'] = c.sha256(ROOT.parent / 'PhysiCelldFBA/config/Ecoli_core.xml')
    assert checks['reference_model_hash'] == checks['current_source_binary_log_hashes']['config/Ecoli_core.xml']
    (HERE / 'validation.json').write_text(json.dumps(checks, indent=2) + '\n')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, metric, title in zip(axes.flat[:2], ['pressure_mean', 'severe_overlap_pairs'],
                                ['Mean reported pressure', 'Pairs exceeding 20% overlap']):
        times = sorted({float(t) for t, _ in old_cells})
        for data, project, label in [(old_cells, 'PhysiCell-X', 'Before'),
                                     (new_cells, 'PhysiCell-X', 'After'),
                                     (new_cells, 'PhysiCelldFBA', 'Reference')]:
            ax.plot(times, [float(data[str(t), project][metric]) for t in times], 'o-', label=label)
        ax.set(xlabel='Time (min)', ylabel=title)
        ax.grid(alpha=.25)
    axes[0, 0].legend()
    for ax, substrate in zip(axes.flat[2:], ['oxygen', 'glucose']):
        for data, label in [(old_fields, 'Before'), (new_fields, 'After')]:
            times = sorted(float(t) for t, s in data if s == substrate and 0 < float(t) <= 180)
            ax.plot(times, [100 * float(data[str(t), substrate]['integrated_change_relative_difference'])
                           for t in times], 'o-', label=label)
        ax.axhline(0, color='black', lw=.8)
        ax.set(xlabel='Time (min)', ylabel=f'{substrate.capitalize()} depletion difference (%)')
        ax.grid(alpha=.25)
    axes[1, 0].legend()
    fig.savefig(HERE / 'before_after.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    xy_extent = max(np.abs(s['labels']['position'][:2]).max() for s in [a, b]) * 1.05
    vmax = max(c.scalar(s, 'pressure').max() for s in [a, b])
    for ax, snap, label in zip(axes, [a, b], ['PhysiCell-X', 'PhysiCelldFBA']):
        pos = snap['labels']['position']
        scatter = ax.scatter(pos[0], pos[1], c=c.scalar(snap, 'pressure'), s=2,
                             vmin=0, vmax=vmax, cmap='magma')
        ax.axvline(0, color='gray', lw=.5)
        ax.set(title=f'{label}, {final_time:g} min', xlabel='x (µm)', ylabel='y (µm)',
               xlim=(-xy_extent, xy_extent), ylim=(-xy_extent, xy_extent), aspect='equal')
    fig.colorbar(scatter, ax=axes, label='Reported pressure (dimensionless)')
    fig.savefig(HERE / 'spatial_pressure_full.png', dpi=180)
    plt.close(fig)
    print('Before/after tables, figures and validation passed.')


if __name__ == '__main__':
    main()
