# AGENTS.md

gmx2qmmm runs QM/MM calculations (single point, optimisation, scan) by coupling
GROMACS (MM) with an external QM code. The repo is mid-refactor (v2): expect
legacy code, thin docs and few tests. Open work is tracked in `PROGRESS.md`.

## Layout

- `src/gmx2qmmm/` — package. `cli.py` (entry point `gmx2qmmm`) → `app.py` (job dispatch)
  - `jobs/` — `singlepoint`, `optimisation`, `mm` (GROMACS), `qm` (QM backends)
  - `generators/` — topology, geometry, point-charge field (PCF), energies
  - `json_files/` — default input parameters, mass map
  - `correction_database/` — correction data shipped as package data
- `examples/{sp,opt}/` — runnable inputs (glycine-serine)
- `test/` — pytest suite (currently one production-style optimisation test)

## QM backends

- Gaussian (`g16`): implemented (`jobs/qm.py`, `QM_gaussian`)
- ORCA: partial
- Turbomole and others: planned
- Planned: a pip-installable QM code (e.g. PySCF) for tests and quickstart, so
  checks run without HPC licenses

## Environment

Target: `pyproject.toml` only, managed with `uv`:

```bash
uv sync --extra test
uv run pytest
```

Not there yet: packaging still lives in `setup.cfg`, and CI uses `env.yaml` (conda)
followed by `pip install -e ".[test]"`. Use that until the migration in `PROGRESS.md` is done.

External programs: GROMACS (`gmx`) is usually available. Gaussian/ORCA are
only available on HPC or specific servers. Never assume a QM code is installed,
and say which checks you could not run.

## Style

Target conventions. Most existing code doesn't follow them yet. Apply them to
code you touch, but don't reformat unrelated code in the same change.

- PEP 8 for naming and layout (enforced with `ruff` once set up, see `PROGRESS.md`)
- Google-style docstrings on public functions, classes and modules
- Type hints on new and changed function signatures
- YAGNI: no options, abstractions or "for later" code nobody needs yet
- Comments explain *why*, not *what*. No chatty or redundant comments, and no
  commented-out code (git keeps the history)
- `pathlib` over `os.path`; `loguru` logger over `print`
- Physical units: state them in names or docstrings (e.g. `energy_hartree`) until pint is introduced

## Rules

- **No behaviour change without a test.** A refactor must keep numerical
  output the same. Compare against reference outputs (e.g.
  `test/optimization_test/ref_output/`) and add a test if none covers the code.
- **Track work in `PROGRESS.md`.** Tick items off when you finish them and add new findings.
- **Ask before deleting** files or code paths, including obvious legacy
  (`*.bkp.py`, committed `env_gmx2qmmm/`, commented-out blocks).
