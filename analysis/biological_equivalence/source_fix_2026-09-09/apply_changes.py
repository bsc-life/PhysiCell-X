"""Apply only the staged changes whose original-file hashes still match."""
import hashlib
import json
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
name = sys.argv[1]
assert name in ['x', 'reference']
manifest = sys.argv[2] if len(sys.argv) > 2 else name + '_changes.json'
assert Path(manifest).name == manifest
records = json.loads((HERE / manifest).read_text())
for record in records:
    path = Path(record['path'])
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None
    assert actual == record['before_sha256'], f'File changed since preparation: {path}'
    assert hashlib.sha256(Path(record['staged']).read_bytes()).hexdigest() == record['after_sha256']
for record in records:
    Path(record['path']).write_bytes(Path(record['staged']).read_bytes())
print('Applied changes:', name, [record['path'] for record in records])
