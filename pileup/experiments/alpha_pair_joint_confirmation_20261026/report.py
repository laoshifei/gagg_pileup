"""Replay the separately frozen confirmation report and notebook delivery."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("validate", "figures", "notebook", "catalog"))
    command = parser.parse_args().command
    folder = Path(__file__).resolve().parent
    freeze = json.loads((folder / "report-freeze.json").read_text())
    for relative, expected in freeze["files"].items():
        actual = hashlib.sha256((folder / relative).read_bytes()).hexdigest()
        assert actual == expected, relative
    assert hashlib.sha256((folder / "freeze.json").read_bytes()).hexdigest() == freeze["inference_freeze_sha256"]
    for name in ("CONFIRM", "JOINT", "NULL", "CONDITIONAL", "RESOLUTION", "PAIR", "COMPONENT"):
        os.environ["GAGG_" + name + "_PROJECT"] = str(folder.parents[2])
    sys.path.insert(0, str(folder / "inference-source"))
    stem = "alpha_pair_joint_confirmation_" + ("notebook" if command in ("notebook", "catalog") else "report")
    spec = importlib.util.spec_from_file_location("pileup.analysis." + stem, folder / "report-source" / (stem + ".py"))
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    getattr(module, "deliver" if command == "notebook" else command)()


if __name__ == "__main__":
    main()
