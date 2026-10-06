"""Assess each completed 30-minute pair, then verify the full source budgets."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_WAIT_POLICY='PASSIVE')
last = -1.


def snapshots(folder):
    result = set()
    for path in folder.glob('output????????.xml'):
        try:
            root = ET.parse(path).getroot()
            names = root.findall('.//simplified_data[@source="PhysiCell"]/filename') + root.findall('microenvironment/domain/data/filename')
            if names and all((folder / node.text).is_file() for node in names):
                result.add(float(root.findtext('metadata/current_time')))
        except ET.ParseError:
            pass  # A newly written XML document may not yet be complete.
    return result


while True:
    folders = [HERE / 'runs' / (name + '_coupled') for name in ['x', 'reference']]
    status = [json.loads((folder / 'run.json').read_text()) for folder in folders]
    assert not any('returncode' in r and not r.get('passed', False) for r in status), 'Simulation failed; inspect run metadata'
    complete = all(r.get('passed', False) for r in status)
    common = snapshots(folders[0]) & snapshots(folders[1])
    endpoint = max(common, default=0.)
    if endpoint >= 30 and (endpoint > last or complete):
        command = [sys.executable, str(ROOT / 'tools/compare_biological_outputs.py'),
                   *map(str, folders), '--out', str(HERE / 'comparison')]
        with (HERE / 'comparison_run.log').open('w') as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([sys.executable, str(HERE / 'summarize_comparison.py')], env=env, check=True)
        last = endpoint
    if complete:
        with (HERE / 'budget_check.log').open('w') as log:
            subprocess.run([sys.executable, str(HERE / 'analyse_budgets.py'), 'coupled'],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        print('Complete: both 180-minute runs, paired snapshots, report and full source budgets verified.', flush=True)
        break
    time.sleep(15)
