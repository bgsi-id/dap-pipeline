#!/usr/bin/env python3
"""List every container image the pipelines reference, for the image-scan matrix.

Sources scanned under pipelines/: `params.<x>_image = '<ref>'` and
`container '<ref>'` in *.nf and *.config, `*_image` keys in params.example.json,
and Dockerfile FROM lines. Null/empty values and references containing `${` are
skipped and reported. A malformed reference fails the run instead of being dropped.
"""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path

REF_RE = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9._-]*(:[0-9]+)?(/[A-Za-z0-9][A-Za-z0-9._-]*)*"
    r"(:[A-Za-z0-9_][A-Za-z0-9_.-]{0,127})?(@sha256:[0-9a-f]{64})?$"
)
DIGEST_RE = re.compile(r"@sha256:[0-9a-f]{64}$")
PARAM_RE = re.compile(r"""(?:params\.)?\b\w*image\s*=\s*(?P<q>['"])(?P<ref>.*?)(?P=q)""")
CONTAINER_RE = re.compile(r"""^\s*container\s+(?P<q>['"])(?P<ref>.*?)(?P=q)""")
NULLISH_RE = re.compile(r"""(?:params\.)?\b\w*image\s*=\s*(null|''|"")\s*$""")
FROM_RE = re.compile(r"^\s*FROM\s+(?:--\S+\s+)*(?P<ref>\S+)(?:\s+AS\s+(?P<alias>\S+))?\s*$", re.IGNORECASE)
GROOVY_SUFFIXES = {".nf", ".config"}


def image_name(ref: str) -> str:
    """Repository path without registry host, tag or digest (the DefectDojo service)."""
    ref = ref.split("@", 1)[0]
    head, sep, tail = ref.rpartition("/")
    if ":" in tail:
        tail = tail.split(":", 1)[0]
    ref = f"{head}{sep}{tail}"
    parts = ref.split("/")
    if len(parts) > 1 and ("." in parts[0] or ":" in parts[0] or parts[0] == "localhost"):
        parts = parts[1:]
    return "/".join(parts)


def _strip_comment(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("//") or stripped.startswith("*") or stripped.startswith("/*"):
        return ""
    return line


def scan_groovy(path: Path) -> list[tuple[str, str]]:
    found = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = _strip_comment(line)
        if not line:
            continue
        if NULLISH_RE.search(line):
            found.append(("", "null"))
            continue
        for regex in (PARAM_RE, CONTAINER_RE):
            match = regex.search(line)
            if match:
                found.append((match.group("ref"), "ref"))
                break
    return found


def scan_json(path: Path) -> list[tuple[str, str]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    found = []
    if isinstance(data, dict):
        for key, value in data.items():
            if key.endswith("_image"):
                found.append((value if isinstance(value, str) else "", "ref"))
    return found


def scan_dockerfile(path: Path) -> list[tuple[str, str]]:
    found = []
    aliases: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        match = FROM_RE.match(line)
        if not match:
            continue
        ref = match.group("ref")
        if ref.lower() == "scratch" or ref in aliases:
            pass
        else:
            found.append((ref, "ref"))
        if match.group("alias"):
            aliases.add(match.group("alias"))
    return found


def discover(root: Path) -> tuple[list[dict], list[str], list[str]]:
    """Return (images, skipped, errors)."""
    images: dict[str, dict] = {}
    skipped: list[str] = []
    errors: list[str] = []
    for path in sorted((root / "pipelines").rglob("*")):
        if not path.is_file():
            continue
        if path.name == "Dockerfile":
            scanner = scan_dockerfile
        elif path.name == "params.example.json":
            scanner = scan_json
        elif path.suffix in GROOVY_SUFFIXES:
            scanner = scan_groovy
        else:
            continue
        rel = path.relative_to(root).as_posix()
        try:
            entries = scanner(path)
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            errors.append(f"{rel}: cannot parse ({exc})")
            continue
        for ref, kind in entries:
            if kind == "null" or not ref.strip():
                skipped.append(f"{rel}: empty or null image (injected at run time)")
                continue
            if "${" in ref or "$" in ref:
                skipped.append(f"{rel}: '{ref}' is resolved at run time")
                continue
            if not REF_RE.match(ref):
                errors.append(f"{rel}: malformed image reference '{ref}'")
                continue
            item = images.setdefault(
                ref,
                {
                    "image_ref": ref,
                    "name": image_name(ref),
                    "pinned": bool(DIGEST_RE.search(ref)),
                    "sources": [],
                },
            )
            if rel not in item["sources"]:
                item["sources"].append(rel)
    return sorted(images.values(), key=lambda i: i["image_ref"]), skipped, errors


def main(argv: list[str]) -> int:
    root = Path(argv[1]) if len(argv) > 1 else Path.cwd()
    images, skipped, errors = discover(root)
    for message in skipped:
        print(f"::notice::skipped {message}")
    for message in errors:
        print(f"::error::{message}")
    unpinned = [i for i in images if not i["pinned"]]
    for item in unpinned:
        print(f"::warning title=Image not pinned to a digest::{item['image_ref']} (used in {', '.join(item['sources'])})")

    payload = json.dumps(images, separators=(",", ":"))
    output = os.environ.get("GITHUB_OUTPUT")
    if output:
        with open(output, "a", encoding="utf-8") as handle:
            handle.write(f"images={payload}\ncount={len(images)}\n")
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as handle:
            handle.write(f"### Container images referenced by the pipelines ({len(images)})\n\n")
            handle.write("| Image | Pinned to digest | Used in |\n| --- | --- | --- |\n")
            for item in images:
                handle.write(f"| `{item['image_ref']}` | {'yes' if item['pinned'] else '**no**'} | {', '.join(item['sources'])} |\n")
    print(payload)

    if errors:
        return 1
    if unpinned and os.environ.get("ENFORCE_DIGEST_PINS", "false").lower() == "true":
        print(f"::error::{len(unpinned)} image(s) are not pinned to a digest and ENFORCE_DIGEST_PINS is true.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
