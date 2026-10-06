"""Exercise deterministic shared-voxel updates in both fresh project builds."""
import json
import os
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
results = []
for name in ['reference', 'x']:
    command = [str(HERE / ('build_' + name) / 'test-cell-sources')]
    env = dict(os.environ, OMP_WAIT_POLICY='PASSIVE', OPENBLAS_NUM_THREADS='1',
               PHYSICELL_SUBSTRATE_BUDGET_PREFIX=str(HERE / ('test_budget_' + name)),
               PHYSICELL_SUBSTRATE_BUDGET_EVERY='10')
    result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=60)
    results.append(dict(project=name, command=command, returncode=result.returncode,
                        stdout=result.stdout, stderr=result.stderr))
target = HERE / 'test_results.json'
if target.exists() and not (HERE / 'test_results_initial.json').exists():
    (HERE / 'test_results_initial.json').write_bytes(target.read_bytes())
target.write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
assert all(r['returncode'] == 0 for r in results)
