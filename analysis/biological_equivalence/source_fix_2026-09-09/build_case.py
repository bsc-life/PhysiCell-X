"""Build fresh isolated colony/test binaries without replacing project objects."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import shutil
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
name = sys.argv[1]
assert name in ['x', 'reference']
incremental = '--incremental' in sys.argv[2:]
recover_tests = '--recover-tests' in sys.argv[2:]
project = ROOT if name == 'x' else ROOT.parent / 'PhysiCelldFBA'
out = HERE / ('build_' + name)
out.mkdir(exist_ok=True)
definitions = dict(re.findall(r'^(\w+)\s*:=\s*(.*)$',
                             (project / 'Makefile').read_text().replace('\\\n', ' '), re.M))


def expand(key):
    return re.sub(r'\$\((\w+)\)', lambda match: expand(match[1]), definitions[key])


objects = expand('ALL_OBJECTS').split()
flags = ['mpic++' if name == 'x' else 'g++', '-std=c++11', '-O3', '-march=native',
         '-fomit-frame-pointer', '-mfpmath=both', '-fopenmp', '-m64', '-DADDON_PHYSIDFBA',
         '-I' + str(project), '-I' + str(project / 'addons/dFBA/ext/libsbml/include'),
         '-I' + str(project / 'addons/dFBA/ext/coin-or/include')]
calls = []


def run(cmd):
    result = subprocess.run([str(p) for p in cmd], cwd=project, text=True, capture_output=True)
    record = dict(command=[str(p) for p in cmd], returncode=result.returncode,
                  stdout=result.stdout, stderr=result.stderr)
    return record


def compile_object(obj):
    stem = Path(obj).stem
    if stem == 'custom':
        source = project / 'custom_modules/custom.cpp'
    else:
        sources = [project / folder / (stem + '.cpp') for folder in
                   ['BioFVM', 'core', 'modules', 'DistPhy', 'addons/dFBA/src']]
        source = next(p for p in sources if p.exists())
    return run(flags + ['-c', source, '-o', out / obj])


if recover_tests:
    assert all((out / obj).exists() for obj in objects)
    calls = json.loads((out / 'compile.json').read_text())
    rebuild_objects = []
elif incremental:
    # Only for the documented final source-helper correction; all other objects
    # come from this experiment's complete, isolated build.
    assert all((out / obj).exists() for obj in objects)
    (out / 'compile_initial.json').write_bytes((out / 'compile.json').read_bytes())
    (out / 'provenance_initial.json').write_bytes((out / 'provenance.json').read_bytes())
    rebuild_objects = ['PhysiCell_cell_container.o']
else:
    rebuild_objects = objects
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for result in pool.map(compile_object, rebuild_objects):
        calls.append(result)
        if result['returncode']:
            print(result['stderr'], flush=True)
(out / 'compile.json').write_text(json.dumps(calls, indent=2) + '\n')
assert all(r['returncode'] == 0 for r in calls), 'Compilation failed; see compile.json'
coin = project / 'addons/dFBA/ext/coin-or/lib'
suffix = '-static.a' if name == 'x' else '.a'
libraries = [coin / ('libClp' + suffix), coin / ('libCoinUtils' + suffix),
             project / 'addons/dFBA/ext/libsbml/lib/libsbml-static.a',
             '-lxml2', '-lbz2', '-lz', '-llapack']
targets = [(project / 'main.cpp', 'bacterial-colony'),
           (HERE / 'test_cell_sources.cpp', 'test-cell-sources')]
for source, target in targets[1:] if recover_tests else targets:
    if not source.exists():
        continue
    # Linking directly on the mounted project filesystem is unnecessarily slow.
    with tempfile.TemporaryDirectory(prefix='physicell-source-link-') as temporary:
        binary = Path(temporary) / target
        result = run(flags + (['-DSOURCE_TEST_MPI'] if name == 'x' else []) +
                     ['-no-pie', source] + [out / obj for obj in objects] + libraries + ['-o', binary])
        if result['returncode'] == 0:
            shutil.copy2(binary, out / target)
    calls.append(result)
    if result['returncode']: print(result['stderr'], flush=True)
(out / 'compile.json').write_text(json.dumps(calls, indent=2) + '\n')
assert all(r['returncode'] == 0 for r in calls), 'Linking failed; see compile.json'
files = [p for folder in ['BioFVM', 'core', 'modules', 'DistPhy', 'addons/dFBA/src', 'custom_modules']
         for p in (project / folder).glob('*') if p.suffix in ['.h', '.cpp']]
files += [project / 'main.cpp', project / 'config/Ecoli_core.xml', out / 'bacterial-colony',
          out / 'test-cell-sources', HERE / 'test_cell_sources.cpp']
(out / 'provenance.json').write_text(json.dumps({str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                              for p in files}, indent=2) + '\n')
print('Fresh build complete:', name, len(objects), 'objects', flush=True)
