"""Create matched configurations, keeping original configs and outputs intact."""
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
cases = []
for name in ['x', 'reference']:
    project = ROOT if name == 'x' else ROOT.parent / 'PhysiCelldFBA'
    for label, duration, interval in [('smoke', .1, .1), ('smoke_final', .1, .1), ('coupled', 180., 30.)]:
        output = HERE / 'runs' / (name + '_' + label)
        output.mkdir(parents=True, exist_ok=True)
        tree = ET.parse(project / 'output/PhysiCell_settings.xml')
        root = tree.getroot()
        root.find('domain/use_2D').text = 'false'
        root.find('parallel/omp_num_threads').text = '2'
        root.find('overall/max_time').text = str(duration)
        root.find('save/full_data/interval').text = str(interval)
        root.find('save/SVG/enable').text = 'false'
        # The MPI cell writer prepends './'; use a relative folder from output/work.
        # Preserve the first smoke configuration as evidence of that I/O limitation.
        root.find('save/folder').text = str(output) if label == 'smoke' else '..'
        for node in root.findall('.//sbml_filename'):
            node.text = str((project / node.text).resolve())
        config = HERE / (name + '_' + label + '.xml')
        tree.write(config, encoding='utf-8', xml_declaration=True)
        cases.append(dict(name=name, label=label, project=str(project), output=str(output),
                          config=str(config), sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
                          mpi_ranks=2 if name == 'x' else 1, omp_threads=2,
                          minutes=duration, snapshot_interval_min=interval,
                          executable=str(HERE / ('build_' + name) / 'bacterial-colony'),
                          budget_prefix=str(output / 'substrate_budget')))
(HERE / 'runs.json').write_text(json.dumps(cases, indent=2) + '\n')
print('Prepared separate smoke tests and matched 180-minute colony runs.')
