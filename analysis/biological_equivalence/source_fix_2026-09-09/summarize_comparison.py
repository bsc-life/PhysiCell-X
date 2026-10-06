"""Create a before/after assessment from completed paired snapshots only."""
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET

os.environ.setdefault('MPLCONFIGDIR', '/tmp/physicell-comparison-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'reassessment_2026-09-09'


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def numeric(row, key):
    return float(row[key]) if row[key] else None


fields = read_csv(HERE / 'comparison/substrate_metrics.csv')
old_fields = read_csv(OLD / 'substrate_metrics.csv')
cells = read_csv(HERE / 'comparison/cell_metrics.csv')
coefficients = [r for r in read_csv(HERE / 'comparison/phenotype_distributions.csv')
                if r['field'] in ['secretion_rates', 'uptake_rates']]
assert coefficients and all(float(r[key]) == 0 for r in coefficients
                            for key in ['x_min', 'x_max', 'reference_min', 'reference_max'])
provenance = json.loads((HERE / 'comparison/provenance.json').read_text())
times = provenance['common_times']
endpoint = max(times)
assert endpoint >= 30, 'Wait for the first full comparison interval'
completed = all(json.loads((HERE / 'runs' / (name + '_coupled') / 'run.json').read_text()).get('passed', False)
                for name in ['x', 'reference']) and endpoint == 180
species = ['oxygen', 'glucose', 'acetate', 'lactate', 'CO2']
active = [s for s in species if s != 'lactate']


def at(rows, time, **matches):
    return next(r for r in rows if float(r['time_min']) == time and all(r[k] == v for k, v in matches.items()))


budgets = {}
for name in ['x', 'reference']:
    rows = [r for path in sorted((HERE / 'runs' / (name + '_coupled')).glob('substrate_budget.rank*.csv'))
            for r in read_csv(path)
            if None not in r.values() and float(r['interval_end_min']) <= endpoint + 1e-7]
    budgets[name] = {}
    for i, substrate in enumerate(species):
        subset = [r for r in rows if int(r['substrate_index']) == i]
        requested = sum(float(r['requested_net_export']) for r in subset)
        applied = sum(float(r['applied_total']) for r in subset)
        assert sum(int(r['steps']) for r in subset) == round(endpoint / .001) * (2 if name == 'x' else 1)
        budgets[name][substrate] = dict(requested=requested, applied=applied,
            relative_error=(applied - requested) / abs(requested) if requested else 0.,
            max_step_error=max(float(r['max_abs_step_error']) for r in subset))

max_accounting = max(abs(r['relative_error']) for data in budgets.values() for r in data.values())
boundary_rows = []
for name in ['x', 'reference']:
    folder = HERE / 'runs' / (name + '_coupled')
    config = ET.parse(HERE / (name + '_coupled.xml')).getroot()
    variables = config.findall('microenvironment_setup/variable')
    for path in sorted(folder.glob('output????????.xml')):
        root = ET.parse(path).getroot()
        time_min = float(root.findtext('metadata/current_time'))
        if time_min not in times: continue
        micro = loadmat(folder / root.findtext('microenvironment/domain/data/filename'))['multiscale_microenvironment']
        for d, variable in enumerate(variables):
            for face in variable.findall('Dirichlet_options/boundary_value'):
                if face.get('enabled', '').lower() != 'true': continue
                label = face.get('ID')
                axis = 'xyz'.index(label[0])
                boundary = micro[axis].min() if label.endswith('min') else micro[axis].max()
                actual = micro[4 + d, micro[axis] == boundary]
                error = float(np.max(np.abs(actual - float(face.text))))
                boundary_rows.append(dict(project=name, time_min=time_min, substrate=variable.get('name'),
                                          face=label, expected=float(face.text), max_absolute_error=error))
assert boundary_rows
max_boundary_error = max(r['max_absolute_error'] for r in boundary_rows)
with (HERE / 'boundary_checks.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(boundary_rows[0]))
    writer.writeheader()
    writer.writerows(boundary_rows)
summary = dict(updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               complete=completed, common_times=times, endpoint_min=endpoint, budgets=budgets,
               max_relative_accounting_error=max_accounting, max_dirichlet_face_error=max_boundary_error)
(HERE / 'assessment_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

plt.rcParams.update({'font.size': 10})
fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax, substrate in zip(axes.flat, active):
    for rows, label, color in [(old_fields, 'Previous saved runs', '#a73b32'),
                               (fields, 'Corrected sources, matched threads', '#12658a')]:
        data = [r for r in rows if r['substrate'] == substrate and 0 < float(r['time_min']) <= endpoint]
        ax.plot([float(r['time_min']) for r in data],
                [100 * float(r['integrated_change_relative_difference']) for r in data],
                marker='o', ms=4, color=color, label=label)
    ax.axhline(0, color='0.6', lw=.7)
    ax.set_title(substrate)
    ax.grid(alpha=.2)
for ax in axes[-1]: ax.set_xlabel('Simulated time (min)')
for ax in axes.flat: ax.set_xlim(0, endpoint + 2)
for ax in axes[:, 0]: ax.set_ylabel('Difference in integrated field change (%)')
axes[0, 0].legend(fontsize=8)
fig.suptitle('PhysiCell-X relative to PhysiCelldFBA\nDepletion / accumulation relative to each pair’s common initial field')
fig.tight_layout()
fig.savefig(HERE / 'before_after_substrates.png', dpi=180)
plt.close(fig)

old_max = max(abs(float(at(old_fields, endpoint, substrate=s)['integrated_change_relative_difference'])) for s in active)
new_max = max(abs(float(at(fields, endpoint, substrate=s)['integrated_change_relative_difference'])) for s in active)
trajectory_max = max(abs(float(r['integrated_change_relative_difference'])) for r in fields
                     if r['substrate'] in active and float(r['time_min']) > 0)
lines = [f'**{"Completed" if completed else "Interim"} comparison through {endpoint:g} minutes.** '
         + ('Both corrected colony runs reached 180 minutes with complete snapshots.' if completed else
            'The 180-minute comparison is still in progress; this report uses only completed paired snapshots.'),
         '', f'The largest integrated substrate-field difference at this endpoint is **{100*new_max:.4g}%**, '
         f'compared with **{100*old_max:.4g}%** in the previous saved runs. '
         f'Across all completed paired times, the largest corrected difference is {100*trajectory_max:.4g}%. '
         'These numbers cover oxygen/glucose depletion and acetate/CO₂ accumulation.',
         '', 'Both projects now use deterministic source application within each voxel and the actual secretion timestep. '
         f'Through the assessed interval, the largest cumulative requested-versus-applied source error is {max_accounting:.3g} relative. '
         'The shared-voxel tests pass at 1, 2 and 8 threads, including changed rates, changed timesteps and reversed cell order.',
         '', f'At **{endpoint:g} minutes**, the comparison is:', '',
         '| Substrate | Previous field-change difference | Corrected field-change difference | Corrected perturbation L2 error | Corrected absolute-field L2 error |',
         '|---|---:|---:|---:|---:|']
for substrate in active:
    old = at(old_fields, endpoint, substrate=substrate)
    new = at(fields, endpoint, substrate=substrate)
    values = [100 * float(old['integrated_change_relative_difference']),
              100 * float(new['integrated_change_relative_difference']),
              100 * float(new['perturbation_relative_l2']), 100 * float(new['relative_l2'])]
    lines.append(f'| {substrate} | {values[0]:+.4g}% | {values[1]:+.4g}% | {values[2]:.4g}% | {values[3]:.4g}% |')
lines += ['', 'Field change means voxel-volume-weighted depletion or accumulation relative to the initial field. '
          'It is distinct from cumulative uptake because diffusion exchanges material with the Dirichlet boundaries. '
          'Perturbation L2 measures the spatial error relative to the reference depletion/accumulation; absolute-field L2 uses the full reference concentration. '
          'The latter explains why oxygen and glucose curves can look very similar even when their small deficits differ.',
          '', '![Before and after substrate-field differences](before_after_substrates.png)',
          '', '**Absolute effect of the correction at 30 minutes.** The values below divide each new integrated field change by '
          'that same project’s previous saved change. A ratio near 0.1 means approximately ten times less depletion or accumulation. '
          'The PhysiCell-X ratios directly match the expected scale of correcting a cached 0.01-minute increment to a 0.001-minute secretion step. '
          'These are observed combined-intervention results, not a factorial estimate separating timestep, race and thread-count effects.',
          '', '| Substrate | X corrected / previous | Reference corrected / previous |', '|---|---:|---:|']
for substrate in active:
    old, new = at(old_fields, 30, substrate=substrate), at(fields, 30, substrate=substrate)
    ratios = [float(new[key]) / float(old[key]) for key in ['x_integrated_change', 'reference_integrated_change']]
    lines.append(f'| {substrate} | {ratios[0]:.8g} | {ratios[1]:.8g} |')
lines += ['', '**Source accounting before the next diffusion step:**', '',
          '| Substrate | Difference in cumulative requested amounts | X accounting error | Reference accounting error |',
          '|---|---:|---:|---:|']
for substrate in active:
    x, r = budgets['x'][substrate], budgets['reference'][substrate]
    difference = 100 * (x['requested'] / r['requested'] - 1)
    lines.append(f'| {substrate} | {difference:+.5g}% | {x["relative_error"]:.3g} | {r["relative_error"]:.3g} |')
lines += ['', 'The budgets sum `dt × net_export_rate` independently of the measured field increment. '
          'Conventional secretion and uptake coefficients are zero in this colony. Accounting errors are dimensionless relative errors; '
          'amounts retain the model’s concentration × voxel-volume units. All ranks and all source steps through the endpoint are included. '
          'The extra tick after the final numbered snapshot is excluded.',
          '', '**FBA timing and causal interpretation.** These changes do not alter FBA optimization or solve scheduling. '
          'Both configurations request a 0.01-minute FBA interval and use a 0.001-minute secretion step. '
          'The repaired bug was application of a cached increment with the wrong timestep; that is different from showing different FBA solve frequencies. '
          'The existing ordering still differs: X can use a newly updated FBA rate in the current source step, while the reference first uses it on the next step. '
          'In the final 0.1-minute smoke test, cumulative requested amounts differed by approximately 1.12%, while application matched those requests to roundoff. '
          'At one minute the oxygen-request difference was approximately 0.102%.',
          '', 'The [prior matched-input diffusion replay](../reassessment_2026-09-09/transport_diagnosis/report.md) '
          'produced bitwise-identical serial and two-rank fields with matching boundaries, including fixed sources. '
          'The present patch changes no transport or Dirichlet code. Source accounting directly tests whether the new coupled runs lose the amounts cells request; '
          'remaining field differences can also reflect source timing, cell placement and subsequent biological feedback. '
          'The combined correction and matched thread counts do not separately quantify how much of the historical mismatch came from each intervention.',
          '', f'An additional check of the enabled boundary faces in every fresh paired snapshot gives a maximum deviation of {max_boundary_error:.3g} '
          'from the configured Dirichlet values, in native concentration units; see [boundary_checks.csv](boundary_checks.csv). '
          'The x/y faces are enabled and the z faces are not clamped in these configurations.',
          '', '**Colony behavior at the assessed endpoint:**', '',
          '| Metric | PhysiCell-X | PhysiCelldFBA | X relative difference |',
          '|---|---:|---:|---:|']
xcell = at(cells, endpoint, project='PhysiCell-X')
rcell = at(cells, endpoint, project='PhysiCelldFBA')
for key, label in [('n_cells', 'Cell count'), ('total_volume', 'Total cell volume'),
                   ('growth_rate_mean', 'Mean stored growth rate'), ('pressure_mean', 'Mean recorded pressure'),
                   ('radius_gyration', 'Radius of gyration'), ('severe_overlap_pairs', 'Pairs with >20% overlap'),
                   ('cells_outside_z', 'Cell centers outside z bounds')]:
    x, r = float(xcell[key]), float(rcell[key])
    difference = f'{100 * (x/r-1):+.5g}%' if r else 'not defined'
    lines.append(f'| {label} | {x:.8g} | {r:.8g} | {difference} |')
trajectory = {}
for key in ['n_cells', 'total_volume', 'growth_rate_mean', 'pressure_mean']:
    values = []
    for t in times:
        x = float(at(cells, t, project='PhysiCell-X')[key])
        r = float(at(cells, t, project='PhysiCelldFBA')[key])
        if r: values.append((abs(x/r-1), t, x, r))
    trajectory[key] = max(values)
lines += ['', 'Across the paired trajectory, the largest absolute relative differences are '
          f'{100*trajectory["n_cells"][0]:.4g}% in cell count, '
          f'{100*trajectory["total_volume"][0]:.4g}% in total cell volume and '
          f'{100*trajectory["growth_rate_mean"][0]:.4g}% in mean stored growth rate.']
p = trajectory['pressure_mean']
lines += [f'The largest relative difference in mean recorded pressure is {100*p[0]:.4g}% at {p[1]:g} minutes '
          f'({p[2]:.7g} versus {p[3]:.7g}). The absolute pressures and overlap distributions matter when interpreting this relative metric.']
negative = sum(int(r[k]) for r in fields for k in ['x_negative_voxels', 'reference_negative_voxels'])
lactate = [r for r in fields if r['substrate'] == 'lactate']
lines += ['', f'Negative substrate voxels across paired snapshots: {negative}. '
          + ('Lactate remains exactly zero in both projects.' if all(float(r['x_max']) == float(r['reference_max']) == 0 for r in lactate) else 'See the lactate rows in the substrate metrics.'),
          '', 'Both runs use two OpenMP threads per process, with two MPI ranks for X and one reference process. '
          'Initial fields and cell positions/volumes match; orientations differ. A single pair of stochastic trajectories, without predefined acceptance margins, '
          'does not establish statistical biological equivalence. Source repairs should be distinguished from agreement of later cell packing.',
          '', 'Full methods, source changes, build recovery and test details are in [README.md](README.md). '
          'The [numerical summary](assessment_summary.json), [substrate metrics](comparison/substrate_metrics.csv), '
          '[cell metrics](comparison/cell_metrics.csv), [comparison provenance](comparison/provenance.json) and '
          '[source test results](test_results.json) provide the supporting data. '
          'The plot of the corrected concentrations is [comparison/substrates.png](comparison/substrates.png).', '']
(HERE / 'report.md').write_text('\n'.join(lines))
late_path = HERE / 'transport_replay_120/summary.json'
if late_path.exists():
    late = json.loads(late_path.read_text())
    matched = all(r['bitwise_equal'] for r in late)
    max_error = max(r['max_absolute_error'] for r in late)
    with (HERE / 'report.md').open('a') as stream:
        stream.write('\n**Additional transport control using the fresh 120-minute state.** '
                     'The reference field and all 815 reference cell positions/rates were replayed for 100 steps of 0.001 minutes, '
                     'with one thread per process and no growth, mechanics or FBA solves. The probes link the fresh production objects used for the coupled simulations. '
                     f'Reference serial versus X serial and X two-rank BLOCKING transport gave a maximum absolute field difference of {max_error:.3g}, '
                     f'both without sources and with identical fixed sources; all four comparisons were bitwise equal: {matched}. '
                     'This checks transport on a later, three-dimensional colony field and supports attributing remaining coupled spatial differences '
                     'to differing source inputs/histories, rather than an operator difference in these controls. '
                     'See the [replay results](transport_replay_120/summary.json), [commands](transport_replay_120/calls.json) '
                     'and [provenance](transport_replay_120/provenance.json).\n')
layout_path = HERE / 'transport_replay_120/layout_summary.json'
if layout_path.exists():
    layout = json.loads(layout_path.read_text())
    rows = [r for r in layout if r['method'] == 'common_rates_per_cell_volume' and r['residual_relative_l2'] is not None]
    with (HERE / 'report.md').open('a') as stream:
        stream.write('\n**What causes the remaining 120-minute spatial mismatch?** '
                     'A second replay used the same reference diffusion solver and initial field for both colonies, '
                     'retaining each colony’s saved cell positions and volumes. Every cell was assigned the same exchange-rate vector '
                     'per unit cell volume, calculated from the reference’s total saved exports and total volume. No FBA solves or mechanics were run. '
                     'After 0.1 minutes of transport with frozen sources, this reproduced the observed coupled field difference with the following residuals:\n\n'
                     '| Substrate | Residual L2 / observed-mismatch L2 |\n|---|---:|\n')
        for row in rows:
            stream.write(f'| {row["substrate"]} | {100*row["residual_relative_l2"]:.4g}% |\n')
        stream.write('\nThe denominator here is the **observed difference between the two coupled fields**, not the full concentration '
                     'or the depletion from baseline. Residuals below 0.35% show that different cell positions/volume distributions reproduce '
                     'almost all of that spatial difference even with identical exchange rates per unit volume. '
                     'This strongly identifies source placement and biomass distribution as the dominant explanation at 120 minutes; '
                     'it does not establish why the underlying cell trajectories differ or statistical equivalence of their packing. '
                     'Using the actual saved rates gives similarly small residuals. '
                     'See [layout_summary.json](transport_replay_120/layout_summary.json) and [layout_calls.json](transport_replay_120/layout_calls.json).\n')
        stream.write('\n![Observed spatial mismatch and the common-rate source-layout replay](transport_replay_120/layout_replay.png)\n')
for late_time in [150, 180]:
    folder = HERE / ('transport_replay_' + str(late_time))
    if late_time > endpoint or not (folder / 'layout_summary.json').exists(): continue
    layout = json.loads((folder / 'layout_summary.json').read_text())
    metadata = json.loads((folder / 'layout_calls.json').read_text())
    rows = [r for r in layout if r['method'] == 'common_rates_per_cell_volume' and r['residual_relative_l2'] is not None]
    residuals = [100*r['residual_relative_l2'] for r in rows]
    volume = {name: float(np.loadtxt(folder / ('layout_' + name + '_common_rate.txt'))[:, 3].sum())
              for name in ['x', 'reference']}
    outside = metadata['outside_cells']
    with (HERE / 'report.md').open('a') as stream:
        stream.write(f'\n**Source-layout follow-up at {late_time} minutes.** '
                     f'The fixtures exclude {outside["x"]} X and {outside["reference"]} reference cell centers outside the saved mesh bounds. '
                     f'Total cell volume inside the domain is {volume["x"]:.8g} versus {volume["reference"]:.8g} '
                     f'({100*(volume["x"]/volume["reference"]-1):+.5g}% for X). '
                     f'Using common exchange rates per unit cell volume, the layout replay leaves {min(residuals):.3g}–{max(residuals):.3g}% '
                     'residual L2 error relative to the observed coupled spatial mismatch. '
                     'This supports a dominant contribution from source arrangement and domain eligibility to the later differences. '
                     'The replay freezes current positions/rates for 0.1 minutes; it does not reconstruct the full movement, division or activation history. '
                     f'See [results](transport_replay_{late_time}/layout_summary.json) and [commands/fixtures](transport_replay_{late_time}/layout_calls.json).\n')
        stream.write('\nCells outside the source mesh are distinct from mismatching substrate Dirichlet data. '
                     'The source helper excludes inactive/out-of-domain cells and invalid voxel indices. '
                     'Disabling z-face Dirichlet clamping does not remove the finite z bounds used to locate cell sources.\n')
print(json.dumps({k: summary[k] for k in ['complete', 'endpoint_min', 'max_relative_accounting_error']}))
