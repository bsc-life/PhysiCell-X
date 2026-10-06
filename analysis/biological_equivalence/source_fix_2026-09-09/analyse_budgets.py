"""Summarize every pre-diffusion source budget from the paired colony runs."""
import csv
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
stage = sys.argv[1] if len(sys.argv) > 1 else 'coupled'
cases = [r for r in json.loads((HERE / 'runs.json').read_text()) if r['label'] == stage]
results = []
for case in cases:
    all_rows = []
    config = ET.parse(case['config']).getroot()
    dt = float(config.findtext('overall/dt_diffusion'))
    expected_steps = round(case['minutes'] / dt)
    paths = sorted(Path(case['output']).glob('substrate_budget.rank*.csv'))
    assert len(paths) == case['mpi_ranks'], (case['name'], paths)
    for path in paths:
        rows = list(csv.DictReader(path.open()))
        assert rows, path
        for row in rows:
            assert all(math.isfinite(float(v)) for v in row.values()), row
            row['rank_file'] = path.name
        # The simulation loop performs one additional diffusion tick after its
        # last numbered snapshot; exclude that partial interval from totals.
        complete = [r for r in rows if float(r['interval_end_min']) <= case['minutes'] + 1e-7]
        for substrate in range(len(config.findall('microenvironment_setup/variable'))):
            series = [r for r in complete if int(r['substrate_index']) == substrate]
            assert sum(int(r['steps']) for r in series) == expected_steps, (path, substrate, 'Incomplete budget')
            previous = 0.
            for row in series:
                assert abs(float(row['interval_start_min']) - previous) < 1e-7, (path, 'Gap or overlap')
                previous = float(row['interval_end_min'])
            assert abs(previous - case['minutes']) < 1e-7, (path, 'Wrong endpoint')
        all_rows += complete
    assert all_rows, 'No complete budget intervals'
    for substrate in sorted({int(r['substrate_index']) for r in all_rows}):
        rows = [r for r in all_rows if int(r['substrate_index']) == substrate]
        requested_net = sum(float(r['requested_net_export']) for r in rows)
        requested_total = sum(float(r['requested_total']) for r in rows)
        applied = sum(float(r['applied_total']) for r in rows)
        max_error = max(float(r['max_abs_step_error']) for r in rows)
        scale = max(abs(requested_net), 1.)
        # This colony has no conventional secretion/uptake terms: its requested
        # source change must equal dt times the net exports, independently of
        # the implicit-source prediction inside the instrumentation.
        assert abs(applied - requested_net) < 1e-9 * scale, (case['name'], substrate, applied, requested_net)
        assert max_error < 1e-8, (case['name'], substrate, max_error)
        results.append(dict(project=case['name'], minutes=case['minutes'], substrate_index=substrate,
                            requested_net_export=requested_net, requested_total=requested_total,
                            applied_total=applied, applied_minus_net=applied-requested_net,
                            relative_accounting_error=(applied-requested_net)/abs(requested_net) if requested_net else None,
                            max_abs_step_error=max_error, rank_steps=sum(int(r['steps']) for r in rows)))
(HERE / (stage + '_budget_summary.json')).write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps(results, indent=2))
