"""Stage reviewable, hash-guarded source changes for both working trees."""
import difflib
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def prepare(project, name):
    replacements = {}
    path = project / 'BioFVM/BioFVM_basic_agent.cpp'
    data = path.read_bytes()
    pattern = rb'\tif\( volume_is_changed \)\r?\n\t\{\r?\n\t\tset_internal_uptake_constants\(dt\);\r?\n\t\tvolume_is_changed = false;\r?\n\t\}'
    newline = b'\r\n' if b'\r\n' in data else b'\n'
    replacement = newline.join([
        b'\t// Rates and the secretion timestep can change independently of volume.',
        b'\t// In particular, dFBA may have cached an intracellular-step increment.',
        b'\tset_internal_uptake_constants(dt);'])
    changed, count = re.subn(pattern, lambda _: replacement, data)
    assert count == 1, (path, count)
    replacements[path] = changed

    path = project / 'core/PhysiCell_cell_container.cpp'
    data = path.read_bytes()
    newline = b'\r\n' if b'\r\n' in data else b'\n'
    include = b'#include "PhysiCell_cell.h"' + newline
    assert data.count(include) == 1
    data = data.replace(include, include + b'#include "PhysiCell_cell_sources.h"' + newline)
    pattern = (rb'\t#pragma omp parallel for[^\r\n]*\r?\n'
               rb'\tfor\( int i=0; i < \(\*all_cells\).size\(\); i\+\+ \)[^\r\n]*\r?\n'
               rb'\t\{\r?\n\t\tif\( \(\*all_cells\)\[i\]->is_out_of_domain == false \)[^\r\n]*\r?\n'
               rb'\t\t\{\r?\n(?:\t\t\tCell\* pC = \(\*all_cells\)\[i\];\r?\n)?'
               rb'\t\t\t(?:pC|\(\*all_cells\)\[i\])->phenotype.secretion.advance\( \(\*all_cells\)\[i\], \(\*all_cells\)\[i\]->phenotype , diffusion_dt_ \);[^\r\n]*\r?\n'
               rb'\t\t\}\r?\n\t\}')
    matches = list(re.finditer(pattern, data))
    assert len(matches) == (2 if name == 'x' else 1), (path, len(matches))
    for index, match in reversed(list(enumerate(matches))):
        rank = b', world.rank' if name == 'x' and index == 1 else b''
        data = data[:match.start()] + b'\tadvance_cell_sources(*all_cells, t, diffusion_dt_' + rank + b');' + data[match.end():]
    replacements[path] = data
    replacements[project / 'core/PhysiCell_cell_sources.h'] = (ROOT / 'core/PhysiCell_cell_sources.h').read_bytes()
    path = project / 'Makefile'
    data = path.read_bytes()
    pattern = rb'(PhysiCell_cell_container\.o: [^\r\n]+)'
    data, count = re.subn(pattern, rb'\1 ./core/PhysiCell_cell_sources.h', data)
    assert count == 1
    replacements[path] = data
    manifest, patch = [], []
    for path, data in replacements.items():
        relative = path.relative_to(project)
        original = path.read_bytes() if path.exists() else None
        staged = HERE / 'staged' / name / relative
        staged.parent.mkdir(parents=True, exist_ok=True)
        staged.write_bytes(data)
        manifest.append(dict(path=str(path), staged=str(staged),
                             before_sha256=hashlib.sha256(original).hexdigest() if original is not None else None,
                             after_sha256=hashlib.sha256(data).hexdigest()))
        patch.extend(difflib.unified_diff((original or b'').decode().splitlines(True), data.decode().splitlines(True),
                                         fromfile='a/' + str(relative), tofile='b/' + str(relative)))
    (HERE / (name + '.patch')).write_text(''.join(patch))
    (HERE / (name + '_changes.json')).write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    prepare(ROOT, 'x')
    prepare(ROOT.parent / 'PhysiCelldFBA', 'reference')
