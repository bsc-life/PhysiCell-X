#!/usr/bin/env python3
"""Compare PhysiCell output folders using XML-labelled MAT data (read-only inputs).

Usage: python3 tools/compare_simulations.py output ../PhysiCell/output --out comparison
Requires numpy, scipy, matplotlib. See compare_simulations.md for metric definitions.
"""
import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
from scipy.io import loadmat
from scipy.spatial import cKDTree


def discover(folder, tolerance, warnings):
    """Use XML times, never filename indices, and exclude pre-initial snapshots."""
    records = []
    for p in sorted(folder.glob('*.xml')):
        if not re.fullmatch(r'(initial|final|output\d+)', p.stem):
            continue
        root = ET.parse(p).getroot()
        time = root.find('metadata/current_time')
        if time is None:
            continue
        if time.get('units', 'min') not in ('min', 'minute', 'minutes'):
            raise ValueError(f'{p}: only simulated times in minutes are supported')
        stamp = root.findtext('metadata/last_modified') or root.findtext('metadata/created')
        epoch = datetime.fromisoformat(stamp.replace('Z', '+00:00')).timestamp() if stamp else p.stat().st_mtime
        records.append(dict(path=p, time=float(time.text), epoch=epoch))
    initial = next((r for r in records if r['path'].stem == 'initial'), None)
    if initial:
        stale = [r for r in records if r['epoch'] < initial['epoch'] - 2]
        for r in stale:
            warnings.append(f"Excluded snapshot older than initial: {r['path']}")
        records = [r for r in records if r not in stale]
    else:
        warnings.append(f'{folder}: initial.xml missing; stale-run detection is limited')
    records.sort(key=lambda r: (r['time'], r['path'].stem != 'initial'))
    unique = []
    for r in records:
        if not np.isfinite(r['time']):
            raise ValueError(f"Invalid time in {r['path']}")
        if unique and abs(r['time'] - unique[-1]['time']) <= tolerance:
            prev = unique[-1]
            if r['path'].stem.startswith('output') and prev['path'].stem.startswith('output'):
                raise ValueError(f"Ambiguous duplicate time: {prev['path']} and {r['path']}")
            # Initial/output-zero and final/numbered aliases are allowed only if
            # their labelled cell state and fields really agree.
            a, b = read_snapshot(prev['path']), read_snapshot(r['path'])
            if not same_state(a, b):
                raise ValueError(f"Conflicting snapshots at the same time: {prev['path']} and {r['path']}")
            continue
        unique.append(r)
    if any(b['epoch'] < a['epoch'] - 2 for a, b in zip(unique, unique[1:])):
        raise ValueError(f'{folder}: timestamps run backwards with simulated time; possible mixed runs')
    if not unique:
        raise ValueError(f'No simulation snapshot XML files found in {folder}')
    return unique


def match_times(a, b, tolerance):
    pairs, unmatched_a, unmatched_b = [], [], []
    i = j = 0
    while i < len(a) and j < len(b):
        delta = a[i]['time'] - b[j]['time']
        if abs(delta) <= tolerance:
            pairs.append((a[i], b[j])); i += 1; j += 1
        elif delta < 0:
            unmatched_a.append(a[i]['time']); i += 1
        else:
            unmatched_b.append(b[j]['time']); j += 1
    return pairs, unmatched_a + [r['time'] for r in a[i:]], unmatched_b + [r['time'] for r in b[j:]]


def read_snapshot(path):
    root = ET.parse(path).getroot()
    nodes = [n for n in root.findall('.//simplified_data') if n.find('labels') is not None]
    domains = root.findall('microenvironment/domain')
    if len(nodes) != 1 or len(domains) != 1:
        raise ValueError(f'{path}: expected one labelled cell population and one microenvironment domain')
    node, domain = nodes[0], domains[0]
    cp, mp = path.parent / node.findtext('filename'), path.parent / domain.findtext('data/filename')
    cells = loadmat(cp)['cells']
    micro = loadmat(mp)['multiscale_microenvironment']
    if not np.isfinite(cells).all() or not np.isfinite(micro).all():
        raise ValueError(f'{path}: non-finite MAT values; fix/inspect outputs before comparison')
    labels, units = {}, {}
    for lab in node.findall('labels/label'):
        name = lab.text.strip()
        start, size = int(lab.get('index')), int(lab.get('size'))
        if name in labels:
            name += '@' + str(start)
        if start < 0 or size < 1 or start + size > cells.shape[0]:
            raise ValueError(f'{path}: invalid label {name}')
        labels[name] = cells[start:start + size]
        units[name] = lab.get('units', '')
    if 'position' not in labels or labels['position'].shape[0] != 3:
        raise ValueError(f'{path}: missing three-component position label')
    if 'ID' in labels and len(np.unique(labels['ID'])) != cells.shape[1]:
        raise ValueError(f'{path}: duplicate cell IDs')
    order = np.lexsort((micro[0], micro[1], micro[2]))
    micro = micro[:, order]
    axes = [np.unique(v) for v in micro[:3]]
    shape = tuple(len(v) for v in axes[::-1])
    if np.prod(shape) != micro.shape[1] or np.any(np.all(np.diff(micro[:3], axis=1) == 0, axis=0)):
        raise ValueError(f'{path}: expected a complete Cartesian mesh with unique voxels')
    if np.any(micro[3] <= 0):
        raise ValueError(f'{path}: non-positive voxel volumes')
    variables = sorted(domain.findall('variables/variable'), key=lambda v: int(v.get('ID')))
    if [int(v.get('ID')) for v in variables] != list(range(len(variables))) or micro.shape[0] != 4 + len(variables):
        raise ValueError(f'{path}: substrate IDs/MAT rows do not agree')
    fields = {v.get('name'): micro[4 + i] for i, v in enumerate(variables)}
    if len(fields) != len(variables):
        raise ValueError(f'{path}: duplicate substrate names')
    return dict(path=path, labels=labels, units=units, n=cells.shape[1], xyz=micro[:3],
                weights=micro[3], axes=axes, shape=shape, fields=fields,
                field_units={v.get('name'): v.get('units', '') for v in variables},
                parameters={v.get('name'): {p.tag: p.text for p in v.findall('physical_parameter_set/*')
                                          if p.tag in ('diffusion_coefficient', 'decay_rate')} for v in variables},
                bounds=np.fromstring(domain.findtext('mesh/bounding_box'), sep=' '),
                files=[path, cp, mp])


def same_state(a, b):
    if a['labels'].keys() != b['labels'].keys() or a['n'] != b['n'] or a['fields'].keys() != b['fields'].keys():
        return False
    if not np.array_equal(a['xyz'], b['xyz']):
        return False
    orders = [np.argsort(scalar(s, 'ID')) if 'ID' in s['labels'] else np.lexsort(s['labels']['position'][::-1])
              for s in (a, b)]
    return (all(np.array_equal(a['labels'][k][:, orders[0]], b['labels'][k][:, orders[1]]) for k in a['labels'])
            and all(np.array_equal(a['fields'][k], b['fields'][k]) for k in a['fields']))


def scalar(s, key):
    return s['labels'].get(key, np.full((1, s['n']), np.nan))[0]


def mean(v):
    return float(np.mean(v)) if len(v) else np.nan


def quantile(v, q):
    return float(np.quantile(v, q)) if len(v) else np.nan


def percent(delta, reference):
    return 100 * delta / abs(reference) if reference else np.nan


def bin_mean(distance, values, edges):
    ids = np.searchsorted(edges, distance, side='right') - 1
    good = (ids >= 0) & (ids < len(edges) - 1) & np.isfinite(values)
    counts = np.bincount(ids[good], minlength=len(edges) - 1)
    sums = np.bincount(ids[good], weights=values[good], minlength=len(edges) - 1)
    return np.divide(sums, counts, out=np.full(len(edges) - 1, np.nan), where=counts > 0)


def nearest_field(s, positions):
    indices = []
    for axis, coord in zip(s['axes'], positions.T):
        right = np.searchsorted(axis, coord).clip(0, len(axis) - 1)
        left = (right - 1).clip(0)
        indices.append(np.where(abs(coord - axis[left]) < abs(coord - axis[right]), left, right))
    return np.ravel_multi_index(tuple(indices[::-1]), s['shape'])


def cell_metrics(s, args):
    pos = s['labels']['position'].T
    radial = np.linalg.norm(pos - args.center, axis=1)
    dead = scalar(s, 'dead') > .5
    phase = scalar(s, 'current_phase')
    live = ~dead
    rate = scalar(s, args.growth_field)
    known_rate = np.isfinite(rate)
    radius = scalar(s, 'radius')
    if len(pos) > 1:
        dist, neighbor = cKDTree(pos).query(pos, k=2)
        nn = dist[:, 1]
        overlap = np.maximum(0, radius + radius[neighbor[:, 1]] - nn)
    else:
        nn = overlap = np.array([])
    metrics = dict(n_cells=s['n'], live_count=int(live.sum()), dead_count=int(dead.sum()),
                   apoptotic_count=int(np.sum(dead & (phase == 100))),
                   necrotic_count=int(np.sum(dead & np.isin(phase, [101, 102, 103]))),
                   growing_count=int(np.sum(live & (rate > args.growth_threshold))) if known_rate.any() else np.nan,
                   arrested_count=int(np.sum(live & (rate <= args.growth_threshold))) if known_rate.any() else np.nan,
                   growth_mean=mean(rate[live]), total_volume=float(scalar(s, 'total_volume').sum()),
                   dead_volume=float(scalar(s, 'total_volume')[dead].sum()),
                   radial_p50=quantile(radial, .5), radial_p95=quantile(radial, .95),
                   radial_max=float(radial.max()) if len(radial) else np.nan,
                   vertical_p95=quantile(abs(pos[:, 2] - args.center[2]), .95),
                   nearest_neighbor_mean=mean(nn), nearest_neighbor_overlap_mean=mean(overlap))
    if 'dead' not in s['labels']:
        for k in ('live_count', 'dead_count', 'apoptotic_count', 'necrotic_count', 'growing_count', 'arrested_count', 'growth_mean', 'dead_volume'):
            metrics[k] = np.nan
    if 'current_phase' not in s['labels']:
        metrics['apoptotic_count'] = metrics['necrotic_count'] = np.nan
    for key in ('pressure', 'total_volume', 'oncoprotein'):
        metrics[key + '_mean'] = mean(scalar(s, key))
        metrics[key + '_p95'] = quantile(scalar(s, key), .95)
    for i, axis in enumerate('xyz'):
        metrics['centroid_' + axis] = mean(pos[:, i])
        metrics['spread_' + axis] = float(np.std(pos[:, i])) if len(pos) else np.nan
    return metrics


def field_comparison(a, b, name, args):
    av, bv = a['fields'][name], b['fields'][name]
    error = av - bv
    weights = b['weights']
    loc = np.argmax(abs(error))
    row = dict(rmse=float(np.sqrt(np.average(error ** 2, weights=weights))),
               relative_l2_percent=percent(np.sqrt(np.sum(error ** 2 * weights)), np.sqrt(np.sum(bv ** 2 * weights))),
               max_absolute_error=float(abs(error[loc])), **{f'max_error_{k}': float(b['xyz'][i, loc]) for i, k in enumerate('xyz')})
    radial = np.linalg.norm(b['xyz'].T - args.center, axis=1)
    for prefix, s, values in [('simulation', a, av), ('reference', b, bv)]:
        for key, val in dict(mean=np.average(values, weights=weights), min=values.min(), max=values.max(),
                             negative_voxels=np.sum(values < 0), core_mean=mean(values[radial < args.core_radius]),
                             cell_exposure_mean=mean(values[nearest_field(s, s['labels']['position'].T)])).items():
            row[prefix + '_' + key] = float(val)
    grad_error = grad_ref = 0.
    for dimension, coords in enumerate(b['axes'][::-1]):
        if len(coords) < 2:
            continue
        ga = np.gradient(av.reshape(a['shape']), coords, axis=dimension)
        gb = np.gradient(bv.reshape(b['shape']), coords, axis=dimension)
        grad_error += np.sum((ga - gb).ravel() ** 2 * weights)
        grad_ref += np.sum(gb.ravel() ** 2 * weights)
    row['gradient_relative_l2_percent'] = percent(np.sqrt(grad_error), np.sqrt(grad_ref))
    return row


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, keys)
        writer.writeheader(); writer.writerows(rows)


def slug(name):
    return re.sub(r'[^a-zA-Z0-9_.-]+', '_', name).strip('_') or 'field'


def plot_lines(plt, path, times, panels, title, labels):
    fig, axes = plt.subplots(int(np.ceil(len(panels) / 3)), 3, figsize=(15, 4 * int(np.ceil(len(panels) / 3))), squeeze=False, layout='constrained')
    for ax, (ylabel, series) in zip(axes.flat, panels):
        for label, vals in zip(labels, series):
            ax.plot(times, vals, '.-', label=label)
        if not any(np.isfinite(v).any() for v in map(np.asarray, series)):
            ax.text(.5, .5, 'Not available in saved outputs', ha='center', transform=ax.transAxes)
        ax.set(xlabel='Simulated time (min)', ylabel=ylabel)
        ax.grid(alpha=.25)
    for ax in list(axes.flat)[len(panels):]:
        ax.set_visible(False)
    axes.flat[0].legend()
    fig.suptitle(title)
    fig.savefig(path, dpi=160); plt.close(fig)


def spatial_plots(plt, out, selected, args, names, edges):
    times = [x[0] for x in selected]
    # These plots deliberately show true thin cell slices, not all cells projected.
    for metric, title in [('pressure', 'Pressure'), ('dead', 'Dead state'), (args.growth_field, 'Growth / cycle rate')]:
        available = [scalar(s, metric) for _, a, b in selected for s in (a, b) if metric in s['labels'] and s['n']]
        if not available:
            continue
        vmax = max(float(v.max()) for v in available); vmin = min(float(v.min()) for v in available)
        fig, axs = plt.subplots(len(selected), 2, figsize=(11, 4 * len(selected)), squeeze=False, layout='constrained')
        for row, (t, a, b) in enumerate(selected):
            for col, s in enumerate((a, b)):
                p = s['labels']['position'].T
                mask = abs(p[:, 2] - args.slice_z) <= args.slice_width / 2
                plot = axs[row, col].scatter(p[mask, 0], p[mask, 1], c=scalar(s, metric)[mask], s=7, cmap='viridis', vmin=vmin, vmax=vmax)
                axs[row, col].axvline(args.interface_x, color='black', ls=':', lw=.8)
                axs[row, col].set(title=f'{args.labels[col]}, t={t:g} min', xlabel='x (µm)', ylabel='y (µm)',
                                  xlim=(args.center[0] - args.spatial_radius, args.center[0] + args.spatial_radius),
                                  ylim=(args.center[1] - args.spatial_radius, args.center[1] + args.spatial_radius), aspect='equal')
        fig.colorbar(plot, ax=axs, label=f"{title} ({selected[-1][1]['units'].get(metric, '')})")
        fig.suptitle(f'{title}: cells within {args.slice_width / 2:g} µm of z={args.slice_z:g}; common color scale')
        fig.savefig(out / f'spatial_cells_{slug(metric)}.png', dpi=160); plt.close(fig)
    for name in names:
        for plane in ('xy', 'xz'):
            fig, axs = plt.subplots(len(selected), 3, figsize=(16, 4 * len(selected)), squeeze=False, layout='constrained')
            slices = []
            for t, a, b in selected:
                axis = 2 if plane == 'xy' else 1
                desired = args.slice_z if plane == 'xy' else args.center[1]
                index = int(np.argmin(abs(a['axes'][axis] - desired)))
                data = [np.take(s['fields'][name].reshape(s['shape']), index, axis=2 - axis) for s in (a, b)]
                slices.append((data, a['axes'][axis][index]))
            low = min(v.min() for data, _ in slices for v in data)
            high = max(v.max() for data, _ in slices for v in data)
            errmax = max(np.max(abs(data[0] - data[1])) for data, _ in slices) or 1e-12
            for row, ((t, a, b), (data, coord)) in enumerate(zip(selected, slices)):
                for col, values in enumerate([*data, data[0] - data[1]]):
                    other = 1 if plane == 'xy' else 2
                    artist = axs[row, col].pcolormesh(a['axes'][0], a['axes'][other], values, shading='nearest',
                                                     cmap='coolwarm' if col == 2 else 'viridis',
                                                     vmin=-errmax if col == 2 else low, vmax=errmax if col == 2 else high)
                    axs[row, col].axvline(args.interface_x, color='black', ls=':', lw=.7)
                    axs[row, col].set(title=f'{(args.labels + ["Simulation − reference"])[col]}, t={t:g}\n{"z" if plane == "xy" else "y"}={coord:g} µm',
                                      xlabel='x (µm)', ylabel=f'{"y" if plane == "xy" else "z"} (µm)',
                                      xlim=(args.center[0] - args.spatial_radius, args.center[0] + args.spatial_radius),
                                      ylim=(args.center[other] - args.spatial_radius, args.center[other] + args.spatial_radius), aspect='equal')
                    if row == 0:
                        fig.colorbar(artist, ax=axs[:, col], label=a['field_units'][name])
            fig.suptitle(name + ': spatial concentration and difference')
            fig.savefig(out / f'spatial_diffusion_{slug(name)}_{plane}.png', dpi=160); plt.close(fig)
        fig, axs = plt.subplots(1, 2, figsize=(14, 5), layout='constrained')
        for j, (t, a, b) in enumerate(selected):
            for col, s in enumerate((a, b)):
                r = np.linalg.norm(s['xyz'].T - args.center, axis=1)
                axs[0].plot((edges[:-1] + edges[1:]) / 2, bin_mean(r, s['fields'][name], edges),
                            ls='-' if col == 0 else '--', color=f'C{j}', label=f'{args.labels[col]}, {t:g} min')
                cylinder = np.linalg.norm(s['xyz'][1:].T - args.center[1:], axis=1) <= args.spatial_radius
                values = s['fields'][name].reshape(s['shape'])
                mask = cylinder.reshape(s['shape'])
                count = mask.sum(axis=(0, 1))
                average = np.divide(np.where(mask, values, 0).sum(axis=(0, 1)), count,
                                    out=np.full(len(s['axes'][0]), np.nan), where=count > 0)
                axs[1].plot(s['axes'][0], average, ls='-' if col == 0 else '--', color=f'C{j}')
        axs[0].set(xlabel='3D radius from specified center (µm)', ylabel=f'{name} ({a["field_units"][name]})')
        axs[0].legend(fontsize=8)
        axs[1].axvline(args.interface_x, color='black', ls=':')
        axs[1].set(xlabel='x (µm)', ylabel='Mean concentration in cylinder', xlim=(args.center[0] - args.spatial_radius, args.center[0] + args.spatial_radius))
        for ax in axs: ax.grid(alpha=.25)
        fig.suptitle(f'{name}: radial and x profiles (cylinder radius {args.spatial_radius:g} µm)')
        fig.savefig(out / f'profiles_diffusion_{slug(name)}.png', dpi=160); plt.close(fig)


def run(args):
    warnings = []
    args.simulation, args.reference, args.out = [p.resolve() for p in (args.simulation, args.reference, args.out)]
    if any(args.out == p or p in args.out.parents for p in (args.simulation, args.reference)):
        raise ValueError('--out must be outside both input folders')
    if args.out.exists() and any(args.out.iterdir()):
        raise ValueError(f'{args.out} is not empty; choose a new --out directory to avoid mixing reports')
    ra, rb = [discover(p, args.time_tolerance, warnings) for p in (args.simulation, args.reference)]
    pairs, ua, ub = match_times(ra, rb, args.time_tolerance)
    if not pairs:
        raise ValueError('No matched simulated times')
    times = [a['time'] for a, b in pairs]
    if abs(times[0]) > args.time_tolerance:
        warnings.append('No matched time-zero snapshot: initial conditions cannot be verified')
    if ua or ub: warnings.append('Unmatched times were omitted; see provenance.json')
    if args.spatial_times:
        select = sorted({min(range(len(times)), key=lambda i: abs(times[i] - t)) for t in args.spatial_times})
    else:
        select = sorted({0, min(1, len(times) - 1), len(times) // 2, len(times) - 1})
    args.out.mkdir(parents=True, exist_ok=True)
    metrics, fields, profiles, manifest, selected = [], [], [], [], []
    initial_info = {}
    for index, (ar, br) in enumerate(pairs):
        a, b = read_snapshot(ar['path']), read_snapshot(br['path'])
        if not np.array_equal(a['xyz'], b['xyz']) or not np.array_equal(a['weights'], b['weights']) or not np.array_equal(a['bounds'], b['bounds']):
            raise ValueError(f'{times[index]} min: meshes differ; interpolation is deliberately not implicit')
        names = sorted(a['fields'].keys() & b['fields'].keys())
        if not names: raise ValueError('No common substrates')
        for name in names:
            if a['field_units'][name] != b['field_units'][name]: raise ValueError(f'{name}: incompatible concentration units')
        for key in a['units'].keys() & b['units'].keys():
            if a['units'][key] != b['units'][key]: raise ValueError(f'{key}: incompatible cell units')
        for s in (a, b):
            if args.growth_field in s['labels'] and s['labels'][args.growth_field].shape[0] != 1:
                raise ValueError(f'{args.growth_field}: growth field must be scalar')
        if index == 0:
            expected_names = names
            if len({slug(k) for k in names}) != len(names):
                raise ValueError('Substrate names collide after conversion to plot filenames')
            initial_info['same_complete_state'] = same_state(a, b)
            initial_info['same_fields'] = all(np.array_equal(a['fields'][k], b['fields'][k]) for k in names)
            if 'ID' in a['labels'] and 'ID' in b['labels']:
                ia, ib = np.argsort(scalar(a, 'ID')), np.argsort(scalar(b, 'ID'))
                initial_info['same_IDs'] = np.array_equal(scalar(a, 'ID')[ia], scalar(b, 'ID')[ib])
                if initial_info['same_IDs']:
                    initial_info['different_labels_by_ID'] = [k for k in a['labels'].keys() & b['labels'].keys()
                        if not np.array_equal(a['labels'][k][:, ia], b['labels'][k][:, ib])]
            r = np.linalg.norm(b['labels']['position'].T - args.center, axis=1)
            if args.spatial_radius is None:
                args.spatial_radius = float(r.max() * 1.2) if len(r) and r.max() else float(max(b['bounds'][3:] - b['bounds'][:3]) / 2)
            if args.core_radius is None: args.core_radius = args.spatial_radius / 3
            edges = np.linspace(0, args.spatial_radius, args.radial_bins + 1)
            for tag, s in zip(args.labels, (a, b)):
                missing = [k for k in ['dead', 'current_phase', 'pressure', 'radius', 'total_volume', args.growth_field] if k not in s['labels']]
                if missing: warnings.append(f'{tag}: missing labels {missing}; associated metrics will be unavailable')
            if a['fields'].keys() != b['fields'].keys(): warnings.append('Only common substrates are compared')
        if names != expected_names: raise ValueError('Substrate schema changes during the run')
        if a['parameters'] != b['parameters']: warnings.append('Diffusion/decay parameters differ in XML metadata')
        for tag, s in zip(('simulation', 'reference'), (a, b)):
            metrics.append(dict(time_min=ar['time'], run=tag, **cell_metrics(s, args)))
            for p in s['files']:
                stat = p.stat()
                manifest.append(dict(path=str(p), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns))
            r = np.linalg.norm(s['labels']['position'].T - args.center, axis=1)
            quantities = {'pressure': scalar(s, 'pressure'), 'dead_fraction': scalar(s, 'dead'),
                          'growth_rate': np.where(scalar(s, 'dead') < .5, scalar(s, args.growth_field), np.nan)}
            for key, values in quantities.items():
                vals = bin_mean(r, values, edges)
                for mid, value in zip((edges[:-1] + edges[1:]) / 2, vals):
                    profiles.append(dict(time_min=ar['time'], run=tag, quantity=key, radius_um=mid, value=value))
        for name in names:
            fields.append(dict(time_min=ar['time'], substrate=name, **field_comparison(a, b, name, args)))
        if index in select: selected.append((ar['time'], a, b))
        print(f'Compared t={ar["time"]:g} min: {a["n"]} / {b["n"]} cells', flush=True)
    # Input mtime checks detect concurrent overwrites, without claiming source/build provenance.
    for entry in manifest:
        stat = Path(entry['path']).stat()
        if stat.st_size != entry['bytes'] or stat.st_mtime_ns != entry['mtime_ns']:
            raise ValueError(f'Input changed during analysis: {entry["path"]}')
    provenance = dict(simulation=str(args.simulation), reference=str(args.reference), matched_times=times,
                      matched_files=[{'simulation': str(a['path']), 'reference': str(b['path']),
                                      'simulation_time': a['time'], 'reference_time': b['time']} for a, b in pairs],
                      simulation_only_times=ua, reference_only_times=ub, first_matched_state=initial_info,
                      spatial_times=[t for t, _, _ in selected], center_um=args.center.tolist(),
                      core_radius_um=args.core_radius, spatial_radius_um=args.spatial_radius,
                      interface_x_um=args.interface_x, growth_field=args.growth_field,
                      growth_threshold=args.growth_threshold, files=manifest, run_metadata={})
    for tag, folder in [('simulation', args.simulation), ('reference', args.reference)]:
        data = {}
        for filename in ('random_seed.txt', 'PhysiCell_settings.xml'):
            p = folder / filename
            if p.exists():
                raw = p.read_bytes()
                data[filename] = dict(sha256=hashlib.sha256(raw).hexdigest())
                if filename.endswith('.txt'): data[filename]['value'] = raw.decode().strip()
            else: warnings.append(f'{tag}: no archived {filename}')
        provenance['run_metadata'][tag] = data
    seeds = [provenance['run_metadata'][tag].get('random_seed.txt', {}).get('value') for tag in ('simulation', 'reference')]
    if seeds[0] != seeds[1]: warnings.append(f'Recorded seeds differ: {seeds[0]} versus {seeds[1]}')
    warnings.append('One run per simulator: descriptive agreement only, not statistical equivalence. Source-code causes are not inferred by this script.')
    provenance['warnings'] = list(dict.fromkeys(warnings))
    write_csv(args.out / 'cell_metrics.csv', metrics)
    write_csv(args.out / 'diffusion_metrics.csv', fields)
    write_csv(args.out / 'radial_cell_profiles.csv', profiles)
    cell_units, field_units = a['units'], a['field_units']
    deltas = []
    for a, b in zip(metrics[::2], metrics[1::2]):
        for key in sorted(a.keys() - {'time_min', 'run'}):
            d = a[key] - b[key]
            deltas.append(dict(time_min=a['time_min'], metric=key, simulation=a[key], reference=b[key],
                               absolute_difference=d, relative_difference_percent=percent(d, b[key])))
    write_csv(args.out / 'metric_differences.csv', deltas)
    (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    with tempfile.TemporaryDirectory(prefix='physicell-mpl-') as cache:
        os.environ.setdefault('MPLCONFIGDIR', cache)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        def panels(spec):
            return [(label, [[r[key] for r in metrics if r['run'] == tag] for tag in ('simulation', 'reference')]) for key, label in spec]
        plot_lines(plt, args.out / 'timeline_mechanics.png', times, panels([
            ('radial_p95', 'Cell-center radius, p95 (µm)'), ('vertical_p95', '|z − center z|, p95 (µm)'),
            ('total_volume', 'Total cell volume (µm³)'), ('pressure_mean', 'Mean saved pressure'),
            ('nearest_neighbor_mean', 'Nearest-neighbor distance (µm)'), ('nearest_neighbor_overlap_mean', 'Nearest-neighbor overlap (µm)')]), 'Mechanics: whole population', args.labels)
        plot_lines(plt, args.out / 'timeline_phenotype.png', times, panels([
            ('n_cells', 'Total cells'), ('live_count', 'Live cells'), ('apoptotic_count', 'Apoptotic cells present'),
            ('necrotic_count', 'Necrotic cells present'), ('growing_count', f'Live cells: {args.growth_field} > threshold'),
            ('arrested_count', 'Live cells at/below growth threshold'), ('growth_mean', f'Mean live {args.growth_field} ({cell_units.get(args.growth_field, "")})'),
            ('dead_volume', 'Total dead-cell volume (µm³)'), ('oncoprotein_mean', 'Mean oncoprotein')]), 'Phenotype: whole population (counts are not cumulative events)', args.labels)
        for name in names:
            rows = [r for r in fields if r['substrate'] == name]
            panel = [(label, [[r[prefix + '_' + key] for r in rows] for prefix in ('simulation', 'reference')])
                     for key, label in [('mean', 'Domain mean'), ('min', 'Domain minimum'), ('core_mean', f'Core mean (r < {args.core_radius:g} µm)'),
                                        ('cell_exposure_mean', 'Concentration at cell positions'), ('negative_voxels', 'Negative voxels'), ('max', 'Domain maximum')]]
            plot_lines(plt, args.out / f'timeline_diffusion_{slug(name)}.png', times, panel, f'{name} ({field_units[name]})', args.labels)
            plot_lines(plt, args.out / f'timeline_errors_{slug(name)}.png', times,
                       [(k.replace('_', ' '), [[r[k] for r in rows]]) for k in ['relative_l2_percent', 'max_absolute_error', 'gradient_relative_l2_percent']],
                       f'{name}: simulation versus reference', ['Error'])
        spatial_plots(plt, args.out, selected, args, names, edges)
        for quantity in ('pressure', 'dead_fraction', 'growth_rate'):
            fig, axs = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
            matrices = [np.array([r['value'] for r in profiles if r['run'] == tag and r['quantity'] == quantity]).reshape(len(times), args.radial_bins)
                        for tag in ('simulation', 'reference')]
            finite = np.concatenate([m[np.isfinite(m)] for m in matrices])
            if not len(finite): plt.close(fig); continue
            for ax, matrix, label in zip(axs, matrices, args.labels):
                art = ax.pcolormesh((edges[:-1] + edges[1:]) / 2, times, matrix, shading='nearest', vmin=finite.min(), vmax=finite.max())
                ax.set(title=label, xlabel='3D radius (µm)', ylabel='Simulated time (min)')
            fig.colorbar(art, ax=axs, label=quantity)
            fig.suptitle('Radial evolution: ' + quantity + ' (empty bins are blank)')
            fig.savefig(args.out / f'radial_evolution_{quantity}.png', dpi=160); plt.close(fig)
    report = ['# Simulation comparison', '', f'Compared {len(times)} matched times, {times[0]:g}–{times[-1]:g} min.', '',
              '| Final metric | Simulation | Reference | Δ | Δ (%) |', '|---|---:|---:|---:|---:|']
    for row in deltas:
        if row['time_min'] == times[-1]:
            report.append(f'| {row["metric"]} | {row["simulation"]:.6g} | {row["reference"]:.6g} | {row["absolute_difference"]:+.6g} | {row["relative_difference_percent"]:+.4g} |')
    report += ['', '## Interpretation limits', ''] + ['- ' + w for w in provenance['warnings']]
    report += ['', 'All population metrics use the full MAT population. Cell maps are thin slices. The dotted x line is a user-selected diagnostic location, not an automatically detected MPI boundary.',
               '', 'Diffusion summaries and physical cell-position exposure can conceal incorrect runtime cell ownership; use the spatial difference maps and x profiles. Actual fluxes, cumulative division/death events, solver state and code provenance cannot be recovered from these snapshots alone.']
    (args.out / 'report.md').write_text('\n'.join(report) + '\n')
    print(f'Graphs, CSV data and report: {args.out}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('simulation', type=Path)
    p.add_argument('reference', type=Path)
    p.add_argument('--out', type=Path, default=Path('comparison'))
    p.add_argument('--labels', nargs=2, default=['Simulation', 'Reference'], metavar=('SIM', 'REF'))
    p.add_argument('--time-tolerance', type=float, default=1e-4)
    p.add_argument('--center', nargs=3, type=float, default=[0, 0, 0], metavar=('X', 'Y', 'Z'))
    p.add_argument('--spatial-radius', type=float, help='Plot radius in µm; default 1.2 × first reference cell radius extent')
    p.add_argument('--core-radius', type=float, help='Core radius in µm; default spatial radius / 3')
    p.add_argument('--radial-bins', type=int, default=30)
    p.add_argument('--slice-z', type=float, default=0)
    p.add_argument('--slice-width', type=float, default=20, help='Thickness of cell z slices in µm')
    p.add_argument('--interface-x', type=float, default=0, help='Diagnostic vertical line in spatial plots (µm)')
    p.add_argument('--spatial-times', nargs='+', type=float, help='Requested times; nearest matched snapshots are used')
    p.add_argument('--growth-field', default='current_cycle_phase_exit_rate', help='XML scalar label; can be growth_rate for metabolic models')
    p.add_argument('--growth-threshold', type=float, default=1e-12)
    args = p.parse_args()
    args.center = np.asarray(args.center)
    if args.radial_bins < 2 or args.time_tolerance < 0 or args.slice_width <= 0 or any(v is not None and v <= 0 for v in (args.spatial_radius, args.core_radius)):
        p.error('Radii/slice width must be positive, radial bins >= 2, time tolerance >= 0')
    try:
        run(args)
    except (ValueError, KeyError, OSError, ET.ParseError) as e:
        p.exit(1, f'Comparison failed: {e}\n')


if __name__ == '__main__':
    main()
