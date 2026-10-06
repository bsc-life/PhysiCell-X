#!/usr/bin/env python3
"""Compile focused production-code probes in /tmp; leave simulation code untouched."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
records = []


def run(args, expected=0):
    result = subprocess.run([str(a) for a in args], cwd=ROOT, text=True, capture_output=True)
    records.append(dict(command=[str(a) for a in args], returncode=result.returncode,
                        stdout=result.stdout, stderr=result.stderr))
    if result.returncode != expected:
        raise RuntimeError(f'Unexpected exit {result.returncode}, expected {expected}: '
                           f'{args}\n{result.stdout}\n{result.stderr}')
    return result.stdout


def main():
    with tempfile.TemporaryDirectory(prefix='physicell-code-probes-') as tmp:
        temp = Path(tmp)
        mesh_source = ROOT / 'BioFVM/BioFVM_mesh.cpp'
        flags = ['mpic++', '-std=c++17', '-O1', '-fopenmp', '-ffunction-sections',
                 '-fdata-sections', '-Wl,--gc-sections', '-IBioFVM']
        mesh_probe = HERE / 'probe_mesh_neighbors.cpp'
        run(flags + [mesh_probe, mesh_source, 'BioFVM/BioFVM_vector.cpp', '-o', temp/'mesh_original'])
        original = run([temp/'mesh_original', '8'], expected=1)
        divisible = run([temp/'mesh_original', '10'])
        # Match the rounding used by resize() in a temporary source copy.
        corrected = mesh_source.read_text()
        for axis, index in [('x', 0), ('y', 1), ('z', 2)]:
            expression = f'(bounding_box[{index+3}]-bounding_box[{index}])/d{axis}'
            old = f'int {axis}_nodes = {expression};'
            new = f'int {axis}_nodes = (int) ceil(1e-16 + {expression});'
            assert corrected.count(old) == 1
            corrected = corrected.replace(old, new)
        (temp/'mesh_corrected.cpp').write_text(corrected)
        run(flags + [mesh_probe, temp/'mesh_corrected.cpp', 'BioFVM/BioFVM_vector.cpp',
                     '-o', temp/'mesh_corrected'])
        repaired = run([temp/'mesh_corrected', '8'])
        # Compile Cell afresh so the filter being tested is the current source.
        cell_flags = ['mpic++', '-std=c++17', '-O1', '-fopenmp', '-DADDON_PHYSIDFBA',
                      '-Iaddons/dFBA/ext/libsbml/include', '-Iaddons/dFBA/ext/coin-or/include']
        run(cell_flags + ['-c', 'core/PhysiCell_cell.cpp', '-o', temp/'cell.o'])
        objects = sorted(p for p in ROOT.glob('*.o') if p.name != 'PhysiCell_cell.o')
        assert objects, 'This filter probe requires the existing bacterial-colony object files.'
        run(cell_flags + [HERE/'probe_mechanics_filter.cpp', temp/'cell.o'] + objects + [
            'addons/dFBA/ext/coin-or/lib/libClp-static.a',
            'addons/dFBA/ext/coin-or/lib/libCoinUtils-static.a',
            'addons/dFBA/ext/libsbml/lib/libsbml-static.a',
            '-lxml2', '-lbz2', '-lz', '-llapack', '-o', temp/'filter'])
        filtered = run([temp/'filter'])
        summary = ('Original mesh, actual 8 µm slab (missing-neighbor failure expected):\n'
                   + original + '\nOriginal mesh, 10 µm divisible control:\n' + divisible
                   + '\nTemporary rounding correction, original 8 µm slab:\n' + repaired
                   + '\nProduction voxel filter, zero versus populated radius cache:\n' + filtered)
        (HERE/'code_probe_results.txt').write_text(summary)
        hashes = {}
        for filename in ['BioFVM/BioFVM_mesh.cpp', 'core/PhysiCell_cell.cpp',
                         'core/PhysiCell_cell_container.cpp', 'core/PhysiCell_standard_models.cpp',
                         'core/PhysiCell_cell_container_mechanics.cpp', 'custom_modules/custom.cpp',
                         'main.cpp']:
            hashes[filename] = hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()
        (HERE/'code_probe_provenance.json').write_text(json.dumps(
            dict(source_sha256=hashes, calls=records), indent=2)+'\n')
        print(summary)


if __name__ == '__main__':
    main()
