"""Check completed budget intervals while the coupled simulations run."""
import csv
import datetime
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
now = datetime.datetime.now(datetime.timezone.utc)
budgets = {}
for name in ['x', 'reference']:
    folder = HERE / 'runs' / (name + '_coupled')
    metadata = json.loads((folder / 'run.json').read_text())
    elapsed = (now - datetime.datetime.fromisoformat(metadata['started_utc'])).total_seconds()
    ranks = []
    data = []
    for path in sorted(folder.glob('substrate_budget.rank*.csv')):
        rows = list(csv.DictReader(path.open()))
        rows = [r for r in rows if None not in r.values()]
        ranks.append(max((float(r['interval_end_min']) for r in rows), default=0.))
        data.extend(rows)
    budgets[name] = data
    errors = [abs(float(r['applied_total']) - float(r['requested_net_export'])) /
              max(abs(float(r['requested_net_export'])), 1.) for r in data]
    print(name, 'wall_minutes=', round(elapsed / 60, 2),
          'simulated_minutes=', [round(t, 3) for t in ranks],
          'max_interval_accounting_error=', max(errors, default=None),
          'finished=', metadata.get('passed', False))

if all(budgets.values()):
    end = min(max(float(r['interval_end_min']) for r in data) for data in budgets.values())
    requested = {name: sum(float(r['requested_net_export']) for r in data
                          if int(r['substrate_index']) == 0 and float(r['interval_end_min']) <= end + 1e-7)
                 for name, data in budgets.items()}
    print('Common oxygen requested difference through', round(end, 3), 'min:',
          (requested['x'] / requested['reference'] - 1) * 100, '%')
