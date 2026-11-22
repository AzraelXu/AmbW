# CLI Output Style Guide

Single reference for terminal output in Bash/Python CLI projects. Business
code must call the shared helpers; do not write ANSI escapes directly.

Modules:

| Language | File | Usage |
|---|---|---|
| Bash | `scripts/lib/ui.sh` | `source "${PROJECT_ROOT}/scripts/lib/ui.sh"` |
| Python | `scripts/sch_ui.py` | `import sch_ui` (add `scripts/` to `sys.path` if needed) |

Both implementations use the same colors, labels, stream routing, and
automatic downgrade rules. They are standard-library only (no colorama/rich):
no network install path, no extra dependency in batch/job environments.

## Color table

| Purpose | SGR | Notes |
|---|---|---|
| FATAL label | `1;97;41` | bold bright white on red |
| FATAL/ERROR body | `1;31` | red bold, whole line |
| WARNING | `1;33` | yellow bold; label only |
| OK / PASS / SUCCESS | `1;32` | green bold; label only |
| INFO / key values / commands / paths | `36` | cyan |
| DEBUG / footnotes / rules | `90` | gray |
| Section headers | `1` | bold, default foreground |
| Reset | `0` | after every colored span |

Rules:

- FATAL and ERROR color the **whole line** red bold. The FATAL label itself is
  white on red.
- WARNING, OK/PASS, INFO, DEBUG, SKIP, FAIL, DRY-RUN color only the leading
  label (and, for Next steps, the command token). The body stays default.
- Plain output, data values, tables, and file contents are never colored.
- Choose colors that read on dark and light terminals; avoid bright blue, dark
  yellow, and bright cyan.

## Modes and precedence

```
--color=auto|always|never   (command line, highest)
SCH_COLOR=auto|always|never (environment default)
auto                        (fallback)
```

`--color` accepts both `--color=always` and `--color always`. `--` stops
option parsing. A missing or invalid value prints
`Invalid --color value: ... (use auto, always, or never)` to stderr and exits 2.

In `auto`, colors are decided **per output stream**:

- stdout helpers check fd 1 is a TTY;
- stderr helpers check fd 2 is a TTY;
- colors are off when the target stream is not a TTY, or `NO_COLOR` is
  non-empty, or `CI` is non-empty, or `TERM=dumb`.

`always` overrides these conditions; `never` disables colors globally.
`cmd > out.log` keeps warnings colored on the terminal; `cmd 2> err.log` writes
no ANSI into `err.log`.

### Bash pre-scan

```bash
source "${PROJECT_ROOT}/scripts/lib/ui.sh"
ui_parse_color_args "$@" || exit $?
set -- "${UI_ARGS[@]+"${UI_ARGS[@]}"}"
```

`ui_parse_color_args` fills `UI_ARGS`, removes `--color`, sets/exports
`SCH_COLOR`, and returns 2 on a bad value. The
`"${UI_ARGS[@]+"${UI_ARGS[@]}"}"` form is required for Bash < 4.4 under
`set -u`.

### Python pre-scan

```python
import sch_ui
args = sch_ui.extract_color_arg(sys.argv[1:])
```

`extract_color_arg` sets the mode, exports `SCH_COLOR`, and exits 2 on a bad
value. Run it before `argparse`.

Aliases/wrappers that exec the script pass `--color` through naturally. Parent
commands export `SCH_COLOR`, so children use the same mode.

## API

Both languages expose:

```
ok  info  warn  error  fatal  debug  header  rule  next_steps
```

Low-level text-preserving emitters:

- Bash: `ui_emit LEVEL TEXT` (stream chosen by level), `ui_emit_out LEVEL TEXT`
  (always stdout)
- Python: `sch_ui.emit(LEVEL, TEXT)`, `sch_ui.emit_out(LEVEL, TEXT)`

`ui_emit` does not add a prefix. It recognizes existing labels:
`[WARN]`, `WARNING :`, `WARNING:`, `WARN:`, `[OK]`, `[PASS]`, `OK:`,
`SUCCESS:`, `[FAIL]`, `FAIL:`, `[INFO]`, `INFO:`, `[DEBUG]`, `DEBUG:`,
`[SKIP]`, `SKIP:`, `[DRY-RUN]`, `[FATAL]`, `FATAL:`, `ERROR:`. If no label
matches, the text is printed unchanged. Use it for legacy prefixes that must
stay byte-identical.

High-level helpers add canonical prefixes: `warn` prints `WARNING: ...`,
`error` prints `ERROR: ...`, `fatal` prints `FATAL: ...`, `ok` prints
`[OK] ...`. `fatal` prints but does **not** exit; callers keep their existing
exit code.

Other helpers:

- `header TEXT`: bold section title.
- `rule [TEXT]`: gray rule; default 60 `=`.
- `debug TEXT`: gray secondary text.
- `ui_kv LABEL VALUE [WIDTH]` / `sch_ui.kv(...)`: preserves the legacy
  `printf '%-Ns: %s'` layout; label bold, value cyan. Python's `kv` has an
  optional `indent` argument.
- `ui_strip TEXT` / `sch_ui.strip(TEXT)`: remove ANSI.
- `next_steps`: see below.

Every Python output function calls `flush=True` and flushes stdout before
writing stderr, preserving order under `cmd 2>&1`.

## Output streams

| Level | Stream |
|---|---|
| FATAL / ERROR / WARNING | stderr |
| OK / INFO / DEBUG / HEADER / RULE / Next steps | stdout |
| `[FAIL]` check results | stdout |
| Report lines that merely contain a level label | unchanged stream (`emit_out`) |

Report lines that are part of a plan/report (a dry-run plan's `[ERROR]` item,
a status report's `WARNING : ETA exceeds ...`) keep their original stream and
field order. Only actual warning/error emissions move to stderr.

If a script writes a log parsed elsewhere, color only the terminal copy. The
line written to `pipeline.log` must be plain, because parsers match
`^\[LEVEL\] timestamp`.

## Next steps

Bash `next_steps STATUS [--footer-line LINE]... STEP...` and Python
`sch_ui.next_steps(status, steps, footer=None)` render:

```
Next steps:
  1) Review generated configs in /work/systems/.
  2) myapp-batch /work/systems/
  3) myapp-status
  4) myapp-results
     (example.conf is a template; directory batch submission skips it
      by default. Pass it explicitly if you really need to run it.)
```

Rules:

1. Numbering is generated from the list; never hardcode.
2. Zero steps means no output at all, including header/footer.
3. `STATUS` is semantic only. The renderer never adds, removes, or filters
   steps; the caller decides whether check/smoke belong in success or failure.
4. Footer lines are gray; first line 5 spaces, continuations 6 spaces.
5. The first token is cyan when it matches `UI_CMD_PREFIXES` (comma-separated,
   default `desmond-`) or looks like `scripts/*.sh` / `*/scripts/*.sh`.
6. Main-flow success: review -> batch -> status -> results. Failure:
   fix/edit -> check -> smoke -> submit -> status -> results.
7. Keep the `example.conf` footnote in the success block.

## Never color these

Even with `--color=always`:

- generated config files and template output;
- managed shell blocks (e.g. `~/.bashrc`), markers, state files;
- log files, job/PBS stdout, JSON reports, result indexes;
- any heredoc, template render, or file write;
- any text another script greps/awks/parses.

Implementation safeguards:

- file writes never go through the UI helpers;
- after redirecting a stage log, call `ui_set_color never` /
  `sch_ui.set_color("never")`;
- captured child processes run with `SCH_COLOR=never`;
- `qsub`/job submission runs with `SCH_COLOR=never`, and generated job scripts
  export `SCH_COLOR=never` after their own redirect;
- `SCH_COLOR` is never frozen into a config.

## Adding a script

1. Source `ui.sh` or import `sch_ui`.
2. Run the `--color` pre-scan before the existing parser.
3. Replace local label printing with the helpers; use `emit` for legacy
   prefixes and `emit_out` for report lines.
4. Keep all `printf > file`, heredocs, and template rendering untouched.
5. If the script writes a log parsed elsewhere, write the plain string and
   color only the terminal copy.
6. Set `UI_CMD_PREFIXES` for the project's command names.
7. Run the skill self-test and the project's full test suite.

## Verification commands

```bash
# No ANSI in a captured log
cmd > out.log 2>&1
LC_ALL=C grep -q $'\033\[' out.log && echo "ANSI leaked" || echo "plain"

# Colors preserved through a pager
cmd --color=always | cat -v | grep -q '\^\[' && echo "colored"

# Automatic downgrade
NO_COLOR=1 cmd | cat -v
CI=true cmd | cat -v
TERM=dumb cmd | cat -v
```

## Python compatibility

The module avoids modern-only syntax and uses `.format()` so it runs on the
older Python bundled with vendor toolchains. It imports only `os`, `re`, and
`sys`.
