#!/usr/bin/env python3
"""Describe paired PhysiCell outputs; no statistical equivalence claim is made.

Requires numpy, scipy and matplotlib. Reads source outputs without modifying them.
Example: python3 tools/compare_biological_outputs.py output ../PhysiCelldFBA/output
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree
from scipy.stats import wasserstein_distance


def snapshots(folder):
    result = {}
    for path in sorted(folder.glob('output????????.xml')):
        root = ET.parse(path).getroot()
        time = float(root.findtext('metadata/current_time'))
        if time in result:
            raise ValueError(f'Duplicate time {time}: {folder}')
        result[time] = path
    return result


def read(path):
    root = ET.parse(path).getroot()
    node = root.find('.//simplified_data[@source="PhysiCell"]')
    cells_path = path.parent / node.findtext('filename')
    micro_path = path.parent / root.findtext('microenvironment/domain/data/filename')
    cells = loadmat(cells_path)['cells']
    micro = loadmat(micro_path)['multiscale_microenvironment']
    assert np.isfinite(cells).all() and np.isfinite(micro).all(), path
    labels, units = {}, {}
    for label in node.findall('labels/label'):
        name = label.text
        # The standard schema repeats elapsed_time_in_phase at indices 8 and 25.
        if name in labels:
            name += '@' + label.attrib['index']
        start, size = int(label.attrib['index']), int(label.attrib['size'])
        labels[name] = cells[start:start + size]
        units[name] = label.attrib.get('units', '')
    ids = labels['ID'][0]
    assert len(np.unique(ids)) == len(ids) and np.equal(ids, ids.astype(int)).all()
    order = np.lexsort((micro[0], micro[1], micro[2]))
    micro = micro[:, order]
    assert len(np.unique(micro[:3].T, axis=0)) == micro.shape[1]
    variables = root.findall('microenvironment/domain/variables/variable')
    names = [v.attrib['name'] for v in variables]
    assert micro.shape[0] == 4 + len(names)
    bounds = np.fromstring(root.findtext('microenvironment/domain/mesh/bounding_box'), sep=' ')
    return dict(path=path, cells=cells, labels=labels, units=units, micro=micro, bounds=bounds,
                names=names, variables=variables, files=[path, cells_path, micro_path])


def scalar(snap, key):
    return snap['labels'][key][0]


def cell_metrics(snap):
    pos = snap['labels']['position'].T
    radii = scalar(snap, 'radius')
    tree = cKDTree(pos)
    pairs = tree.query_pairs(2 * radii.max(), output_type='ndarray')
    dist = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)
    overlap = np.maximum(0, 1 - dist / (radii[pairs[:, 0]] + radii[pairs[:, 1]]))
    positive = overlap > 0
    overlaps = overlap[positive]
    severe = pairs[overlap > .2]
    cross_midplane = pos[severe[:, 0], 0] * pos[severe[:, 1], 0] < 0
    # Independently reconstruct the standard cell-cell pressure from geometry.
    geometric_pressure = np.zeros(len(pos))
    for side in (0, 1):
        np.add.at(geometric_pressure, pairs[:, side], overlap**2 / 0.027288820670331)
    center = pos.mean(axis=0)
    radius_gyration = np.sqrt(np.mean(np.sum((pos - center)**2, axis=1)))
    metrics = dict(n_cells=len(pos), total_volume=scalar(snap, 'total_volume').sum(),
                   dead_count=int(scalar(snap, 'dead').sum()),
                   radius_gyration=radius_gyration,
                   nearest_neighbor_mean=tree.query(pos, k=2)[0][:, 1].mean(),
                   overlap_pairs=int(positive.sum()),
                   overlap_mean=overlaps.mean() if len(overlaps) else 0.,
                   overlap_max=overlaps.max() if len(overlaps) else 0.,
                   severe_overlap_cells=len(np.unique(severe)),
                   severe_overlap_pairs=len(severe),
                   severe_overlap_pairs_crossing_x_zero=int(cross_midplane.sum()),
                   geometric_pressure_mean=geometric_pressure.mean(),
                   geometric_pressure_max=geometric_pressure.max(),
                   z_min=pos[:, 2].min(), z_max=pos[:, 2].max(),
                   cells_outside_z=int(((pos[:, 2] < snap['bounds'][2]) |
                                       (pos[:, 2] > snap['bounds'][5])).sum()),
                   velocity_max=np.linalg.norm(snap['labels']['velocity'], axis=0).max())
    for name in ['pressure', 'total_volume', 'radius', 'growth_rate']:
        values = scalar(snap, name)
        for suffix, value in zip(['mean', 'min', 'p50', 'p95', 'max'],
                                 [values.mean(), values.min(), *np.quantile(values, [.5, .95]), values.max()]):
            metrics[name + '_' + suffix] = value
    for name in ['dead', 'current_phase', 'cycle_model', 'cell_type', 'is_motile']:
        vals, counts = np.unique(scalar(snap, name), return_counts=True)
        metrics[name + '_counts'] = json.dumps(dict(zip(vals.astype(int).astype(str), counts.tolist())))
    for graph in ['cell_neighbor', 'attached_cells', 'spring_attached_cells']:
        path = snap['path'].with_name(snap['path'].stem + '_' + graph + '_graph.txt')
        lines = path.read_text().splitlines() if path.exists() else []
        metrics[graph + '_graph_rows'] = len(lines)
    return metrics


def ratio(numerator, denominator):
    return float(numerator / denominator) if denominator != 0 else None


def field_metrics(a, b, time, baseline):
    assert a['names'] == b['names']
    assert np.array_equal(a['micro'][:4], b['micro'][:4]), 'Mesh mismatch'
    xyz = b['micro'][:3]
    axes = [np.unique(row) for row in xyz]
    shape = tuple(len(axis) for axis in axes[::-1])
    assert np.prod(shape) == xyz.shape[1]
    roi = xyz[0]**2 + xyz[1]**2 <= 25**2
    weights = b['micro'][3]
    rows = []
    for i, name in enumerate(b['names']):
        av, bv = a['micro'][4 + i], b['micro'][4 + i]
        error = av - bv
        da, db = av - baseline[i], bv - baseline[i]
        rmse = np.sqrt(np.mean(error**2))
        grad_a = np.array(np.gradient(av.reshape(shape), *axes[::-1], edge_order=2))
        grad_b = np.array(np.gradient(bv.reshape(shape), *axes[::-1], edge_order=2))
        row = dict(time_min=time, substrate=name, x_mean=av.mean(), reference_mean=bv.mean(),
                   x_min=av.min(), reference_min=bv.min(), x_max=av.max(), reference_max=bv.max(),
                   rmse=rmse, max_absolute_error=np.max(np.abs(error)),
                   relative_l2=ratio(np.linalg.norm(error), np.linalg.norm(bv)),
                   perturbation_relative_l2=ratio(np.linalg.norm(error), np.linalg.norm(db)),
                   gradient_relative_l2=ratio(np.linalg.norm(grad_a-grad_b), np.linalg.norm(grad_b)),
                   x_integrated_change=np.sum(da*weights), reference_integrated_change=np.sum(db*weights),
                   integrated_change_relative_difference=ratio(np.sum(da*weights), np.sum(db*weights)),
                   x_roi_mean=av[roi].mean(), reference_roi_mean=bv[roi].mean(),
                   roi_rmse=np.sqrt(np.mean(error[roi]**2)),
                   x_negative_voxels=int((av < 0).sum()), reference_negative_voxels=int((bv < 0).sum()))
        if row['integrated_change_relative_difference'] is not None:
            row['integrated_change_relative_difference'] -= 1
        rows.append(row)
    return rows


def csv_write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('--out', type=Path, default=Path('analysis/biological_equivalence'))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    files_a, files_b = snapshots(args.project), snapshots(args.reference)
    times = sorted(files_a.keys() & files_b.keys())
    if not times or times[0] != 0:
        raise ValueError('Comparison requires common snapshots including time zero')
    cells, fields, distributions, manifest = [], [], [], []
    for folder in [args.project, args.reference]:
        for name in ['PhysiCell_settings.xml', 'random_seed.txt']:
            path = folder / name
            if path.exists():
                manifest.append(dict(path=str(path.resolve()), bytes=path.stat().st_size,
                                     sha256=sha256(path)))
    baseline = None
    initial = None
    for time in times:
        a, b = read(files_a[time]), read(files_b[time])
        assert a['labels'].keys() == b['labels'].keys(), 'Cell schema mismatch'
        assert a['units'] == b['units'], 'Cell units mismatch'
        assert np.array_equal(a['bounds'], b['bounds']), 'Domain mismatch'
        for va, vb in zip(a['variables'], b['variables']):
            assert va.attrib == vb.attrib, 'Substrate names/units mismatch'
            for parameter in ['diffusion_coefficient', 'decay_rate']:
                tag = 'physical_parameter_set/' + parameter
                assert float(va.findtext(tag)) == float(vb.findtext(tag)), parameter
        if baseline is None:
            assert np.array_equal(a['micro'], b['micro']), 'Initial fields differ'
            baseline = b['micro'][4:].copy()
            ia, ib = np.argsort(scalar(a, 'ID')), np.argsort(scalar(b, 'ID'))
            assert np.array_equal(scalar(a, 'ID')[ia], scalar(b, 'ID')[ib])
            initial = dict(differing_labels=[key for key in a['labels']
                           if not np.array_equal(a['labels'][key][:, ia], b['labels'][key][:, ib])],
                           initial_count=len(ia), identical_initial_fields=True)
        for label, snap in [('PhysiCell-X', a), ('PhysiCelldFBA', b)]:
            cells.append(dict(time_min=time, project=label, **cell_metrics(snap)))
            for path in snap['files']:
                manifest.append(dict(path=str(path.resolve()), bytes=path.stat().st_size,
                                     sha256=sha256(path)))
            for graph in ['cell_neighbor', 'attached_cells', 'spring_attached_cells']:
                path = snap['path'].with_name(snap['path'].stem + '_' + graph + '_graph.txt')
                if path.exists():
                    manifest.append(dict(path=str(path.resolve()), bytes=path.stat().st_size,
                                         sha256=sha256(path)))
        fields.extend(field_metrics(a, b, time, baseline))
        for name in a['labels']:
            if name in ['ID', 'position', 'orientation', 'attack_target']:
                continue
            assert a['labels'][name].shape[0] == b['labels'][name].shape[0]
            for component, (av, bv) in enumerate(zip(a['labels'][name], b['labels'][name])):
                distributions.append(dict(time_min=time, field=name, component=component,
                    units=a['units'][name], x_mean=av.mean(), reference_mean=bv.mean(),
                    relative_mean_difference=ratio(av.mean()-bv.mean(), abs(bv.mean())),
                    wasserstein=wasserstein_distance(av, bv),
                    x_min=av.min(), reference_min=bv.min(), x_max=av.max(), reference_max=bv.max()))
        if time == 30:
            ia, ib = np.argsort(scalar(a, 'ID')), np.argsort(scalar(b, 'ID'))
            if np.array_equal(scalar(a, 'ID')[ia], scalar(b, 'ID')[ib]):
                displacement = np.linalg.norm(a['labels']['position'][:, ia]-b['labels']['position'][:, ib], axis=0)
                initial['position_difference_at_30_min'] = dict(mean=displacement.mean(), max=displacement.max())
        last_a, last_b = a, b
    csv_write(args.out / 'cell_metrics.csv', cells)
    csv_write(args.out / 'substrate_metrics.csv', fields)
    csv_write(args.out / 'phenotype_distributions.csv', distributions)
    (args.out / 'provenance.json').write_text(json.dumps(dict(
        project=str(args.project.resolve()), reference=str(args.reference.resolve()),
        common_times=times, project_only_times=sorted(files_a.keys()-files_b.keys()),
        reference_only_times=sorted(files_b.keys()-files_a.keys()), initial=initial,
        files=manifest), indent=2) + '\n')
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/physicell-comparison-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    metrics = [('n_cells', 'Cell count'), ('total_volume', 'Total cell volume (µm³)'),
               ('pressure_mean', 'Mean reported pressure (dimensionless)'),
               ('radius_gyration', 'Radius of gyration (µm)'),
               ('overlap_max', 'Maximum fractional overlap'), ('growth_rate_mean', 'Mean growth rate (h⁻¹)')]
    for ax, (key, label) in zip(axs.flat, metrics):
        for project in ['PhysiCell-X', 'PhysiCelldFBA']:
            rows = [r for r in cells if r['project'] == project]
            ax.plot(times, [r[key] for r in rows], 'o-', label=project)
        ax.set(xlabel='Time (min)', ylabel=label)
        ax.grid(alpha=.25)
    axs[0, 0].legend()
    fig.savefig(args.out / 'comparison.png', dpi=180)
    plt.close(fig)
    fig, axs = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for ax, name in zip(axs.flat, last_b['names']):
        rows = [r for r in fields if r['substrate'] == name]
        for prefix, label in [('x', 'PhysiCell-X'), ('reference', 'PhysiCelldFBA')]:
            ax.plot(times, [r[prefix+'_roi_mean'] for r in rows], 'o-', label=label)
        ax.set(xlabel='Time (min)', ylabel=name + ' (saved concentration units)')
        ax.grid(alpha=.25)
    axs[0, 0].legend()
    axs[1, 2].axis('off')
    fig.suptitle('Mean concentration within 25 µm of colony origin, all z layers')
    fig.savefig(args.out / 'substrates.png', dpi=180)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, snap, label in zip(axs, [last_a, last_b], ['PhysiCell-X', 'PhysiCelldFBA']):
        pos = snap['labels']['position']
        plot = ax.scatter(pos[0], pos[1], c=scalar(snap, 'pressure'), s=7, vmin=0, vmax=.2, cmap='magma')
        ax.set(title=f'{label}, {times[-1]:g} min', xlabel='x (µm)', ylabel='y (µm)',
               xlim=(-22, 22), ylim=(-22, 22), aspect='equal')
    fig.colorbar(plot, ax=axs, label='Reported pressure (values above 0.2 saturated)')
    fig.savefig(args.out / 'spatial_pressure.png', dpi=180)
    plt.close(fig)
    print(json.dumps(dict(common_times=times, initial=initial, final_cells=cells[-2:],
                          final_fields=fields[-5:]), indent=2))


if __name__ == '__main__':
    main()
