#!/usr/bin/env python3
"""Validate catalog.yaml against the contract DAP enforces when it syncs pipelines.

Fails (exit 1) when an entry is malformed, ids repeat, a path escapes
`pipelines/<dir>`, a descriptor is missing or inconsistent, or an entrypoint
does not exist. Pipeline directories that are not in the catalog are only listed.
"""
from __future__ import annotations

import sys
from pathlib import Path, PurePosixPath

import yaml

REQUIRED_KEYS = ("id", "path", "descriptor")


def _load_yaml(path: Path):
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _safe_pipeline_dir(raw: object) -> str | None:
    """Return an error message when `raw` is not a safe `pipelines/<dir>` path."""
    if not isinstance(raw, str) or not raw:
        return "path must be a non-empty string"
    if raw.startswith(("/", "\\")):
        return f"path '{raw}' must be relative"
    parts = raw.split("/")
    if any(part in ("", ".", "..") for part in parts):
        return f"path '{raw}' must not contain '.', '..' or empty segments"
    if len(parts) != 2 or parts[0] != "pipelines":
        return f"path '{raw}' must be a directory directly under pipelines/"
    return None


def validate(root: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    infos: list[str] = []
    catalog_file = root / "catalog.yaml"
    if not catalog_file.is_file():
        return ["catalog.yaml not found"], infos
    try:
        catalog = _load_yaml(catalog_file)
    except yaml.YAMLError as exc:
        return [f"catalog.yaml is not valid YAML: {exc}"], infos
    if not isinstance(catalog, dict):
        return ["catalog.yaml must be a mapping"], infos
    if not catalog.get("schema_version"):
        errors.append("catalog.yaml: schema_version is missing")
    entries = catalog.get("pipelines")
    if not isinstance(entries, list):
        return errors + ["catalog.yaml: 'pipelines' must be a list"], infos

    seen: set[str] = set()
    listed_dirs: set[str] = set()
    for index, entry in enumerate(entries):
        label = f"pipelines[{index}]"
        if not isinstance(entry, dict):
            errors.append(f"{label}: entry must be a mapping")
            continue
        entry_id = entry.get("id")
        if isinstance(entry_id, str) and entry_id:
            label = f"'{entry_id}'"
        missing = [key for key in REQUIRED_KEYS if not isinstance(entry.get(key), str) or not entry.get(key)]
        if missing:
            errors.append(f"{label}: missing or empty {', '.join(missing)}")
        if not isinstance(entry.get("active"), bool):
            errors.append(f"{label}: 'active' must be a boolean")
        if missing:
            continue
        if entry_id in seen:
            errors.append(f"{label}: duplicate id")
        seen.add(entry_id)

        path_error = _safe_pipeline_dir(entry["path"])
        if path_error:
            errors.append(f"{label}: {path_error}")
            continue
        pipeline_dir = root / entry["path"]
        listed_dirs.add(PurePosixPath(entry["path"]).name)
        if not pipeline_dir.is_dir():
            errors.append(f"{label}: directory {entry['path']} does not exist")
            continue

        descriptor_rel = PurePosixPath(entry["descriptor"])
        if descriptor_rel.is_absolute() or ".." in descriptor_rel.parts:
            errors.append(f"{label}: descriptor '{entry['descriptor']}' must be a relative path without '..'")
            continue
        try:
            descriptor_rel.relative_to(entry["path"])
        except ValueError:
            errors.append(f"{label}: descriptor '{entry['descriptor']}' is not inside {entry['path']}")
            continue
        descriptor_file = root / descriptor_rel
        if not descriptor_file.is_file():
            errors.append(f"{label}: descriptor {entry['descriptor']} does not exist")
            continue
        try:
            descriptor = _load_yaml(descriptor_file)
        except yaml.YAMLError as exc:
            errors.append(f"{label}: descriptor {entry['descriptor']} is not valid YAML: {exc}")
            continue
        if not isinstance(descriptor, dict):
            errors.append(f"{label}: descriptor {entry['descriptor']} must be a mapping")
            continue
        if descriptor.get("id") != entry_id:
            errors.append(
                f"{label}: descriptor id '{descriptor.get('id')}' differs from catalog id '{entry_id}'"
            )
        entrypoint = descriptor.get("entrypoint")
        if not isinstance(entrypoint, str) or not entrypoint:
            errors.append(f"{label}: descriptor has no entrypoint")
        else:
            entry_rel = PurePosixPath(entrypoint)
            if entry_rel.is_absolute() or ".." in entry_rel.parts:
                errors.append(f"{label}: entrypoint '{entrypoint}' must be relative without '..'")
            elif not (pipeline_dir / entrypoint).is_file():
                errors.append(f"{label}: entrypoint {entry['path']}/{entrypoint} does not exist")

    pipelines_dir = root / "pipelines"
    if pipelines_dir.is_dir():
        for child in sorted(pipelines_dir.iterdir()):
            if child.is_dir() and child.name not in listed_dirs:
                infos.append(f"pipelines/{child.name} is not in catalog.yaml (unpublished)")
    return errors, infos


def main(argv: list[str]) -> int:
    root = Path(argv[1]) if len(argv) > 1 else Path.cwd()
    errors, infos = validate(root)
    for info in infos:
        print(f"::notice::{info}")
    for error in errors:
        print(f"::error::{error}")
    if errors:
        print(f"Catalog validation failed with {len(errors)} error(s).")
        return 1
    print("Catalog validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
