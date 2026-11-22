#!/usr/bin/env python3
"""Strict {{KEY}} template renderer for Amber input files.

No shell evaluation, no expression language: every placeholder is replaced by
a literal string supplied with --set or --json. Missing and unused keys are
errors by default, which catches typos before tleap/pmemd does.

Usage:
  python3 render_template.py TEMPLATE --set KEY=VALUE ... [--json FILE]
                              [--out FILE] [--force] [--lenient]
                              [--allow-unused] [--list-keys]

Exit codes: 0 success, 2 usage/validation error.
"""
import argparse
import json
import os
import re
import sys

PLACEHOLDER = re.compile(r"\{\{([A-Za-z_][A-Za-z0-9_]*)\}\}")


def parse_set(items):
    values = {}
    for item in items or []:
        if "=" not in item:
            sys.exit(f"--set expects KEY=VALUE, got: {item}")
        key, value = item.split("=", 1)
        key = key.strip()
        if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", key):
            sys.exit(f"invalid placeholder name: {key!r}")
        values[key] = value
    return values


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("template")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    ap.add_argument("--json", metavar="FILE", help="JSON object of key/value pairs")
    ap.add_argument("--out", metavar="FILE", help="write here instead of stdout")
    ap.add_argument("--force", action="store_true", help="overwrite --out")
    ap.add_argument("--lenient", action="store_true",
                    help="leave unresolved placeholders in place")
    ap.add_argument("--allow-unused", action="store_true",
                    help="do not fail on values the template does not use")
    ap.add_argument("--list-keys", action="store_true",
                    help="print the template's placeholders and exit")
    args = ap.parse_args()

    if not os.path.isfile(args.template):
        sys.exit(f"template not found: {args.template}")
    text = open(args.template, encoding="utf-8").read()
    keys = PLACEHOLDER.findall(text)

    if args.list_keys:
        for key in sorted(set(keys)):
            print(key)
        return

    values = {}
    if args.json:
        with open(args.json, encoding="utf-8") as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            sys.exit(f"--json must contain an object, got {type(data).__name__}")
        values.update({str(k): str(v) for k, v in data.items()})
    values.update(parse_set(args.set))

    needed = set(keys)
    missing = sorted(needed - set(values))
    if missing:
        sys.exit("missing values for: " + ", ".join(missing))

    unused = sorted(set(values) - needed)
    if unused and not args.allow_unused:
        sys.exit("values not used by template: " + ", ".join(unused))

    def repl(match):
        key = match.group(1)
        return values[key]

    rendered = PLACEHOLDER.sub(repl, text)

    if "{{" in rendered and not args.lenient:
        # Catch malformed or unknown placeholder syntax.
        bad = sorted(set(re.findall(r"\{\{[^}]*\}?\}?", rendered)))
        sys.exit("unresolved/malformed placeholders: " + ", ".join(bad))

    if args.out:
        if os.path.exists(args.out) and not args.force:
            sys.exit(f"refusing to overwrite {args.out} (use --force)")
        parent = os.path.dirname(os.path.abspath(args.out))
        os.makedirs(parent, exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(rendered)
    else:
        sys.stdout.write(rendered)


if __name__ == "__main__":
    main()
