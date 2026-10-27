"""Replay reports using their own snapshot and the frozen inference dependencies."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('validate', 'figures', 'notebook'))
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    manifest = json.loads((folder / 'report-freeze.json').read_text())
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest(folder / 'freeze.json') == manifest['inference_freeze_sha256']
    for relative, expected in manifest['files'].items():
        assert digest(folder / relative) == expected, relative
    for name in ('SEARCH', 'JOINT', 'CONFIRM', 'NULL', 'CONDITIONAL', 'RESOLUTION', 'PAIR', 'COMPONENT'):
        os.environ['GAGG_' + name + '_PROJECT'] = str(folder.parents[2])
    sys.path.insert(0, str(folder / 'inference-source'))
    suffix = 'notebook' if args.command == 'notebook' else 'report'
    name = 'pileup.analysis.alpha_pair_weak_search_' + suffix
    spec = importlib.util.spec_from_file_location(name, folder / 'report-source' / (name.rsplit('.', 1)[1] + '.py'))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    getattr(module, 'deliver' if args.command == 'notebook' else args.command)()


if __name__ == '__main__':
    main()
