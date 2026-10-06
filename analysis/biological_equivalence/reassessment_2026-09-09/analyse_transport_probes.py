"""Summarize the transport isolation experiments; raw replay arrays live in /tmp."""
import csv
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MPLCONFIGDIR', '/tmp/physicell-comparison-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / 'transport_diagnosis'
ROOT = HERE.parents[2]


def read(paths):
    values = np.concatenate([np.fromfile(p).reshape(-1, 8) for p in paths])
    assert len(np.unique(values[:, :3], axis=0)) == 250000
    assert np.isfinite(values).all()
    return values[np.lexsort((values[:, 0], values[:, 1], values[:, 2]))]


def main():
    source_summary = []
    for path in sorted(OUT.glob('source_*.csv')):
        source_rows = list(csv.DictReader(path.open()))
        for name in ['oxygen', 'glucose', 'acetate', 'CO2']:
            errors = np.array([float(r['relative_error']) for r in source_rows
                               if r['substrate'] == name])
            source_summary.append(dict(case=path.stem, substrate=name,
                                       mean_error_percent=100 * float(errors.mean()),
                                       minimum_error_percent=100 * float(errors.min()),
                                       maximum_error_percent=100 * float(errors.max())))
    (OUT / 'source_summary.json').write_text(json.dumps(source_summary, indent=2) + '\n')
    results = []
    for replay in [0, 1]:
        reference_path = f'/tmp/biological-transport-reference-replay{replay}_0.bin'
        reference = read([reference_path])
        for solver in ['serial', 'blocking2']:
            if solver == 'serial':
                paths = [f'/tmp/biological-transport-x-replay{replay}_0.bin']
            else:
                prefix = '/tmp/biological-transport-blocking2' + ('-replay1' if replay else '')
                paths = [prefix + f'_{rank}.bin' for rank in [0, 1]]
            values = read(paths)
            assert np.array_equal(values[:, :3], reference[:, :3])
            error = values[:, 3:] - reference[:, 3:]
            results.append(dict(solver=solver, fixed_sources=bool(replay), steps=100,
                                step_minutes=.001, bitwise_identical=np.array_equal(values, reference),
                                max_absolute_error_by_substrate=np.max(abs(error), axis=0).tolist(),
                                paths=paths, reference=reference_path,
                                file_sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                                             for p in paths + [reference_path]}))
    (OUT / 'diffusion_summary.json').write_text(json.dumps(results, indent=2) + '\n')
    rows = list(csv.DictReader((HERE / 'substrate_metrics.csv').open()))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    for name in ['oxygen', 'glucose', 'acetate', 'CO2']:
        data = [r for r in rows if r['substrate'] == name and float(r['time_min']) > 0]
        times = [float(r['time_min']) for r in data]
        axes[0].plot(times, [100 * (float(r['x_roi_mean']) / float(r['reference_roi_mean']) - 1)
                            for r in data], 'o-', label=name)
        axes[1].plot(times, [100 * float(r['integrated_change_relative_difference'])
                            for r in data], 'o-', label=name)
    for ax in axes:
        ax.axhline(0, color='black', lw=.8)
        ax.set(xlabel='Time (min)', ylabel='Difference from reference (%)')
        ax.grid(alpha=.25)
    axes[0].set_title('Concentration: same 25 µm ROI as substrates.png')
    axes[1].set_title('Whole-domain depletion / accumulation from t = 0')
    axes[0].legend()
    fig.savefig(OUT / 'substrate_relative_differences.png', dpi=180)
    plt.close(fig)
    filenames = ['addons/dFBA/src/dfba_intracellular.cpp', 'addons/dFBA/src/dfba_intracellular.h',
                 'addons/dFBA/src/dfba_Model.cpp', 'core/PhysiCell_cell_container.cpp',
                 'core/PhysiCell_phenotype.cpp', 'BioFVM/BioFVM_basic_agent.cpp',
                 'BioFVM/BioFVM_microenvironment.cpp', 'BioFVM/BioFVM_solvers.cpp',
                 'output/PhysiCell_settings.xml', 'config/Ecoli_core.xml', 'main.cpp']
    provenance = {}
    for project in [ROOT, ROOT.parent / 'PhysiCelldFBA']:
        for filename in filenames:
            path = project / filename
            provenance[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        for filename in ['BioFVM_vector.o', 'BioFVM_mesh.o', 'BioFVM_microenvironment.o',
                         'BioFVM_solvers.o', 'BioFVM_agent_container.o', 'BioFVM_utilities.o',
                         'BioFVM_matlab.o', 'pugixml.o']:
            path = project / filename
            provenance[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    for filename in ['transport_probe.cpp', 'cache_timestep_probe.cpp', 'build_transport_probe.py',
                     'analyse_transport_probes.py', 'reproduce_transport_probes.py']:
        path = HERE / filename
        provenance[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    provenance['/tmp/biological-transport-initial.bin'] = hashlib.sha256(
        Path('/tmp/biological-transport-initial.bin').read_bytes()).hexdigest()
    (OUT / 'source_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
