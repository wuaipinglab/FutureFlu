#!/usr/bin/env python3
"""
Shared helpers for the release step layer.

These helpers keep step loading and temporary work paths consistent.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


def load_project_lib(project_root: Path):
    # Load the package-local workflow core for this release tree.

    lib_path = project_root / "scripts" / "workflow_core.py"
    if not lib_path.exists():
        raise FileNotFoundError(f"missing project library: {lib_path}")

    spec = importlib.util.spec_from_file_location("workflow_core", lib_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load project library from {lib_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def component_part_dir(lib, run_root: Path) -> Path:
    # Component slices stay under output-root/_step_outputs, not in published trees.

    del lib  # package-local helper keeps the same call signature for all steps
    path = run_root.parent / "_step_outputs" / run_root.name / "component_parts"
    path.mkdir(parents=True, exist_ok=True)
    return path


def run_root_for_theta(lib, theta: float) -> Path:
    # Match FutureFlu package layout under package outputs/: no theta_* nesting.
    # Nest by theta only for isolated experiment output roots.
    out = Path(lib.OUTPUT_ROOT)
    default_outputs = (lib.ROOT / "outputs").resolve()
    if out.resolve() == default_outputs:
        return out
    return out / lib.theta_label(theta)
