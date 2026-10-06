"""Run the requested matched smoke/coupled cases and record exact commands."""
import concurrent.futures
import datetime
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
stage = sys.argv[1]
assert stage in ['smoke', 'smoke_final', 'coupled']
cases = [r for r in json.loads((HERE / 'runs.json').read_text()) if r['label'] == stage]


def run(case):
    output = Path(case['output'])
    assert not (output / 'final.xml').exists(), 'Refusing to overwrite a completed run'
    work = output / 'work'
    work.mkdir(exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               OMP_WAIT_POLICY='PASSIVE', PHYSICELL_SUBSTRATE_BUDGET_PREFIX=case['budget_prefix'],
               PHYSICELL_SUBSTRATE_BUDGET_EVERY='100' if stage.startswith('smoke') else '1000')
    cmd = [case['executable'], case['config']]
    if case['name'] == 'x':
        cmd = ['mpirun', '--allow-run-as-root', '--oversubscribe', '--bind-to', 'none', '-n', '2'] + cmd
    metadata = dict(case=case, command=cmd, cwd=str(work),
                    started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    environment={k: env[k] for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                      'MKL_NUM_THREADS', 'OMP_WAIT_POLICY', 'PHYSICELL_SUBSTRATE_BUDGET_PREFIX',
                      'PHYSICELL_SUBSTRATE_BUDGET_EVERY']})
    (output / 'run.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Starting', case['name'], stage, case['minutes'], 'simulated minutes', flush=True)
    start = time.monotonic()
    with (output / 'run.log').open('w') as log:
        process = subprocess.Popen(cmd, cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT)
        metadata['pid'] = process.pid
        (output / 'run.json').write_text(json.dumps(metadata, indent=2) + '\n')
        metadata['returncode'] = process.wait()
    metadata['elapsed_seconds'] = time.monotonic() - start
    metadata['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    final = output / 'final.xml'
    metadata['final_time_min'] = float(ET.parse(final).getroot().findtext('metadata/current_time')) if final.exists() else None
    snapshots = sorted(output.glob('output????????.xml'))
    missing = []
    for snapshot in snapshots + ([final] if final.exists() else []):
        root = ET.parse(snapshot).getroot()
        for node in root.findall('.//simplified_data[@source="PhysiCell"]/filename') + root.findall('microenvironment/domain/data/filename'):
            if not (output / node.text).is_file():
                missing.append(str(output / node.text))
    metadata['missing_snapshot_files'] = missing
    metadata['passed'] = metadata['returncode'] == 0 and metadata['final_time_min'] is not None and metadata['final_time_min'] >= case['minutes'] - 1e-6 and bool(snapshots) and not missing
    (output / 'run.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('Finished', case['name'], 'passed=', metadata['passed'], 'wall_seconds=', metadata['elapsed_seconds'], flush=True)
    return metadata


with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(run, cases))
assert all(r['passed'] for r in results), 'A simulation failed; inspect run.json and run.log'
