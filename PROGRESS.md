Roadmap
-------

Rough order. Each phase makes the next one safe. IDs refer to the findings below.

0. **PR #26 (``optimization`` -> ``master``)**: settle the issues the PR introduced or touches
   before merging, so they don't land on ``master``. See "PR #26" below.
1. **Safety net**
   - Fix the CI trigger (H1)
   - Make tests run from any CWD, write only to ``tmp_path``, and stop being circular (T1–T4)
   - Clean up tracked junk and ``.gitignore`` (H2–H4, H6; ask before deleting)
   - Move to a ``pyproject.toml``-only install that works with both ``uv`` and conda; pick one Python version (H5)
   - Lint with ``ruff`` in CI (fix H7). Add formatting later, file by file
2. **Regression baseline**: numeric golden tests for SP energies/forces, PCF and
   generated topology, plus unit tests for pure helpers. Freeze current
   behaviour before touching it. (T5–T7)
3. **Fix the remaining verified bugs**, one per change, each with a test. (B1–B16 not handled in phase 0)
   If a fix changes reference numbers, regenerate the golden files and record why in the commit.
4. **Config**: YAML + Pydantic schema with defaults and validation. It replaces
   ``params.txt`` parsing (``input_handler.py``) and the unused JSON defaults (C1–C5, B1–B3).
   It also settles path handling: one ``Path`` type, resolved against the config
   file or ``work_dir``, never the CWD, including GROMACS executable/share dirs (B4, B6).
5. **QM backend interface**: one base class for writing input, running and parsing.
   Port Gaussian, add PySCF (pip-installable, runs everywhere, used for tests and
   quickstart), then ORCA. Also covers subprocess handling and moving work out of constructors. (A1–A5, B5, B7)
6. **Topology/system rewrite**: replace the hand-rolled ``.top`` parsing with
   [MDParser](https://github.com/janjoswig/MDParser) (``GromacsTopologyParser``: resolves includes and
   ``#ifdef``, typed entries, can write ``.top`` back). Check first:
   - support for ``[ constraints ]``, ``[ virtual_sites* ]``, ``[ cmap ]``, and whether it can write the merged QM/MM ``.top``
   - packaging: PyPI only has 0.0.2, CI is on Travis; it needs a release
   Then: units via pint (U1–U3), fix performance (P1–P8, S4).
7. **Cleanup and docs**: remove dead code (ask first), make the style follow AGENTS.md. (D1–D6)

PR #26 (``optimization``)
-------------------------

Author: AlinaJansen. Status (2026-11-13): open, review "changes requested" (janjoswig), with a Copilot review.
Reconciled against all review comments. *(new)* = found in the audit and not yet raised on the PR.

**Fix on the PR** (small, introduced or touched by the PR; author to fix):

- [ ] H2 Uncommit ``env_gmx2qmmm/``; gitignore ``.DS_Store`` and remove ``examples/opt/.DS_Store``
- [ ] H3 Delete ``optimization_production_test.bkp.py``, ``params.bkp.txt`` and ``oenergy.bkp.txt``. *(new)*: ``test_files/oenergy.txt``/``oforces.txt`` are run outputs that get rewritten on every test run
- [ ] B9 ``databasecorrection`` (47fe0bd): UnboundLocalError when the method *is* in the DB, a redundant second query, ``fetchall()[0]`` IndexError with a misplaced length check, and a temp file that is never deleted
- [ ] B16 Rejected steps: history entries are appended before evaluation and never removed (the PR dropped the ``.pop()`` calls); output is written under the decremented step number
- [ ] B16 *(new)* ``optimisation.py:87`` compares against the raw string ``f_thresh`` (``float()`` is used at 167/180 but not here)
- [ ] B8 ``conjgrad`` calls an undefined function (``optimisation.py:134``); propagator case is mixed, so make it case-insensitive. *(new)*: also raise a clear error for CONJGRAD/BFGS, since ``geometry.py:169/214`` are broken too
- [ ] B6 ``optimisation.py:212`` cleanup uses bare filenames, which fails when ``work_dir`` != CWD. *(new)*: same for QM files (``qm.py:141ff``)
- [ ] U1 ``geometry.py:162`` extra ``* 0.52917721``: explain it with a comment and check consistency (``geometry.py:214`` uses ``0.052917721``)
- [ ] ``mm.py:53`` explain or revert disabling the large-box enlargement; drop the duplicate ``initial``/``large`` attributes; keep box vectors as an array
- [ ] ``mm.py:66`` ``edr_xvg_file`` should probably be ``edrname.with_suffix(".xvg")``; ``*name`` attributes are really paths
- [ ] Small cleanups: debug ``print`` at ``optimisation.py:201`` (delete), stray ``pass`` (``:144``), f-string for ``pcf_filename`` (``:106``)
- [ ] B4 ``check_gmx_paths``: missing separator in ``gmxpath + gmxcmd``. Reviewer asks for ``Path`` objects and clear names: ``gmx_exe`` (absolute, from ``shutil.which``), drop ``gmxpath`` or make it ``gmx_bin_dir``, ``gmxtop_path`` -> ``gmx_top_dir``/``gmx_shared_dir``
- [ ] *(new)* H1 CI trigger is still ``branches: main``
- [ ] *(new)* CI hasn't run since 2026-05-28 (last run passed). d9177bd and 9a2a9d6 have no checks; maybe workflow approval is needed for the pushing account
- [ ] *(new)* d9177bd (water fix: ``-DFLEXIBLE`` instead of SETTLE) is untested and may change MM output drastically for systems with water. Needs a with/without-water comparison before merging

**Propose on a separate branch, merge into the PR** (larger, maintainer to do):

- [ ] T1–T4 Production test rework, per the discussion in the PR conversation: drop the fully mocked test (it passes by construction), add a real reference run plus unit tests for the parts (link correction, displacement, PCF update). Also: CWD-independent, ``tmp_path``, local ``amberGS.ff`` lookup (currently fails outside CI), ``zip`` -> length-checked numeric compare
- [ ] Defer to roadmap phase 4: ``gmx_*`` path renaming beyond the B4 fix, if the author only does the minimal fix
- [ ] Defer to roadmap phase 5: SCAN branching inside ``Optimisation`` (``:138``). A scan should be its own job that wraps optimisations

Miscellaneous items to address
------------------------------

Things not covered by the roadmap.

- [ ] README is outdated: Python 3.6, GROMACS 5.0.2, Gaussian only, ``python gmx2qmmm.py`` usage, input-file table of CLI flags that no longer exist
- [ ] README advertises the SCAN job type, but ``app.py`` doesn't dispatch it. Drop it from the README or implement it
- [ ] No ``LICENSE`` file, although ``setup.cfg`` declares MIT
- [ ] User docs: decide where they live (``gmx2qmmm.github.io/gmx2qmmm_io``?) and how they're built
- [ ] Stale branches: local ``{new_branch_name}``, remotes ``NA_MA_Integration``, ``migration_*``, ``prep_release_v1``, ``refactor/pcf``, ``dev``. Check with owners and prune
- [ ] Version is hardcoded ``2.0.0`` in ``__init__.py``; decide on release and versioning scheme

Findings from codebase audit (2026-11-13)
-----------------------------------------

Line numbers are as of commit d9177bd.

### Hygiene / packaging / CI

- [ ] H1 CI triggers on push to ``main``, but the default branch is ``master``, so only PRs run CI (``.github/workflows/python-app.yml``)
- [ ] H2 Committed virtualenv ``env_gmx2qmmm/`` (points to ``/homes/chemie/jansen``)
- [ ] H3 Tracked junk: ``examples/*/.DS_Store``, ``examples/sp/gro_info.txt``, ``*.bkp.py``/``*.bkp.txt`` under ``test/optimization_test`` (contain ``/afs/...`` paths)
- [ ] H4 ``.gitignore`` ignores ``examples/`` and ``params.txt``; the re-include rules can't work; ``!test/unit/test_files/*`` points to a dir that doesn't exist
- [ ] H5 Metadata split across ``setup.cfg``/``setup.py``/``pyproject.toml``; Python 3.9 (``setup.cfg``) vs 3.11 (``env.yaml``)
- [ ] H6 ``examples/*/params.txt`` hardcode ``/trinity/...`` GROMACS paths
- [ ] H7 Invalid escape sequences in regex strings (``"\S+"``), which emit ``SyntaxWarning`` on Python >= 3.12

### Tests

- [ ] T1 ``test/optimization_test`` fails locally: ``amberGS.ff`` is looked up only in the GROMACS share dir, not next to the topology
- [ ] T2 Tests use paths relative to the CWD (must run from the repo root) and write into the repo (``logfile.log*``, ``oenergy.txt``, ``test_write_g96``). Use ``tmp_path``
- [ ] T3 The optimisation test is circular: the link correction is mocked as ``qm+mm-ref``, so the total matches the reference by construction
- [ ] T4 File comparisons use ``zip()`` (a truncated output still passes) and exact string match (fragile against float formatting)
- [ ] T5 Golden numeric tests on existing ``test/production/ref_output`` (PCF, ``.qmmm.top``, ``.ndx``, energies, forces)
- [ ] T6 Unit tests for pure helpers: ``_helper`` (``_flatten``, ``filter_xyzq``, link-atom coords), g96 round trip, ``propagate_dispvec`` STEEP
- [ ] T7 No coverage for multi-link systems, methods present in the correction DB, ``.gro`` input, or ``-wd`` != CWD

### Verified bugs

- [ ] B1 ``input_handler.py:84ff`` ``split('=')[1]`` truncates values that contain ``=`` (e.g. ``iop(3/33=1)``); a line without ``=`` raises IndexError
- [ ] B2 ``input_handler.py:146`` list values are split but never stored (they stay raw strings)
- [ ] B3 ``input_handler.py:192`` ``ValueError`` is logged and the function returns ``None``, so ``app.py:55`` crashes later
- [ ] B4 ``app.py:177`` auto-detected ``gmxpath`` has no trailing slash, but ``mm.py`` builds ``gmxpath + gmxcmd``, giving a broken command
- [ ] B5 Subprocess return codes ignored (``mm.py``, ``qm.py``); with ``-backup no``, stale files get parsed
- [ ] B6 Mixed CWD vs ``work_dir`` paths (``qm.py:141ff``, ``optimisation.py:212``, ``mm.py:46/111``, ``system.py:96/572``) break any run with ``-wd`` != CWD
- [ ] B7 ``singlepoint.py:67`` any ``qmcommand`` other than ``g16`` leaves the QM job undefined (ORCA crashes with AttributeError)
- [ ] B8 ``optimisation.py:131`` ``conjgrad`` calls a nonexistent method; ``bfgs`` reads an unset attribute; ``geometry.py:169/214`` CONJGRAD/BFGS are wrong (``[] *= 0.0529``). Only STEEP works
- [ ] B9 ``energies.py:138`` correction DB lookup raises UnboundLocalError whenever the method *is* in the DB. Only works today because BP86/STO-3G falls back to CAM-B3LYP silently
- [ ] B10 ``energies.py:394`` link-atom energy correction is applied once after the loop (last link only); the force version applies it per link
- [ ] B11 ``energies.py:353`` vs ``:583`` energy and force loops disagree on the number of linkcorr charges (``len(m2)`` vs ``len(m2)*2``)
- [ ] B12 ``energies.py:441`` ``get_m2charges()`` returns ``None``, then gets indexed (TypeError)
- [ ] B13 ``qm.py:67`` ``range(1, N)`` drops the last atom (off-by-one)
- [ ] B14 ``qm.py:137`` an existing ``.gjf.log`` silently skips the QM step (stale results)
- [ ] B15 ``mm.py:143`` ``cutoff-scheme = group`` was removed in GROMACS 2020; ``mm.py:268`` energy term "11" is hardcoded on stdin
- [ ] B16 ``optimisation.py:87`` compares against unconverted ``f_thresh`` (TypeError for non-integer values); a rejected step is not rolled back (suspected)

### Suspected bugs (multi-link PCF; needs a test system to confirm)

- [ ] S1 ``pcf.py:557`` ragged ``reshape`` fails on numpy >= 1.24 when M1 atoms have different numbers of M2 neighbours
- [ ] S2 ``pcf.py:559-573`` M2 shift charges are written in atom order while the field is grouped by M1, so charges land on the wrong atoms with more than one M1
- [ ] S3 ``pcf.py:168-185`` a ``continue`` skips ``count += 1``, so the M1 lists go out of sync
- [ ] S4 ``system.py:861ff``/``topology.py:1179`` link lists are indexed with flattened Q1 indices, which misaligns if an M1 has more than one QM neighbour
- [ ] S5 ``pcf.py:294`` convergence check compares displacement, not change, so it never fires

### Architecture

- [ ] A1 QM backend chosen by scattered string checks on ``qmcommand`` and ``program``; ``QM`` base class contains Gaussian code; ``QM_gaussian`` doesn't call ``super()``
- [ ] A2 Constructors do all the work (``App``, ``SystemInfo``, ``GenerateTopology``, ``Singlepoint``, ``Optimisation``)
- [ ] A3 Shared mutable state: ``input_dict`` mutated in jobs; step counter lives on ``system``
- [ ] A4 Filename scheme ``jobname[.step].ext`` rebuilt in 4 places; ``qm.py:95`` redirects ``sys.stdout`` to write files
- [ ] A5 ``exit(1)``/``sys.exit`` instead of exceptions throughout; ``logger.remove()`` touches global state

### Config

- [ ] C1 ``json_files/input_parameters_default.json`` is never loaded, and its key names differ from the ones the code reads
- [ ] C2 No defaults and no validation (``Asserter`` is a no-op); a missing key raises KeyError deep inside a job
- [ ] C3 Floats are never converted; only lines that start with ``#`` count as comments
- [ ] C4 ``#include``/``!INCLUDE`` paths resolve against the CWD, not the including file
- [ ] C5 ``base_dir`` (``app.py:39``) is unused and wrong once pip-installed

### Parsing / topology

- [ ] P1 ``.top`` section loops end only on a blank line: the next header is swallowed or raises ``int('[')``
- [ ] P2 ``[ constraints ]``, ``[ virtual_sites* ]``, ``[ cmap ]``, ``[ position_restraints ]`` are silently dropped
- [ ] P3 Section-header regexes are inconsistent (``[dihedrals]`` skipped; ``[moleculetype]`` gives length 0 and shifts offsets)
- [ ] P4 Molecule-name match has no word boundary (``PRO`` matches ``PROT``)
- [ ] P5 ``#ifdef`` emulation can't nest; the block is pasted about 8 times in ``make_large_top`` (~840-line function)
- [ ] P6 Charge regexes fail on ``+0.5`` or scientific notation
- [ ] P7 ``read_gmx_box_vectors`` handles rectangular boxes only; it raises UnboundLocalError without a BOX
- [ ] P8 Performance: the molecule file is re-parsed per molecule instance; ``x in np.array(list)`` in nested loops; ``mm.py:52`` builds an all-pairs N×N×3 distance matrix (out of memory on solvated systems)

### Units

- [ ] U1 ``0.52917721`` hardcoded about 60 times; ``geometry.py:214`` uses ``0.052917721``
- [ ] U2 Unexplained constants: ``0.7409471631``, ``0.5282272551``, ``0.71290813568205`` (two of them are the same value written twice), PCF factors ``0.1944``/``0.0233``
- [ ] U3 ``mm.py:228/307`` kJ/mol to Hartree factors unnamed; units of ``stepsize``/``f_thresh`` undocumented

### Dead code / duplication (ask before deleting)

- [ ] D1 Unused modules: ``charge_shift.py``, ``prepare_pcf_for_shift.py``, ``sum_pcf_tm.py`` (``pcf.py`` is a class copy of them)
- [ ] D2 ``pcf_from_top.py``: only ``checkformol`` is used; ``readg96``/``readgeo`` duplicate ``geometry.py``
- [ ] D3 ``correction_database/database.csv`` is never read (identical to the sqlite)
- [ ] D4 Unused: ``logging_utils.Logger``, ``_helper.create_dict``/``stepper``/``make_xyzq``, ``STEPLIMIT``, ``nmaflag``, many imports
- [ ] D5 "Find moleculetype by name" implemented 5 times; ``.gro`` parsing 3 times; z/v link lookup 4 times in ``energies.py``
- [ ] D6 Large commented-out blocks (``geometry.py:129-213``, ``topology.py:1279``, ``system.py:1023``); ``__main__`` inside a class body (``topology.py:1296``)

Summary of changes related to the migration to v2 (historic)
------------------------------------------------------------

  * Switched to a `setup.cfg`/`pyproject.toml` based installation mechanism
  * Switched to a `src/` based package layout
  * Exposed main CLI command as a proper entry point
    * Separated `gmx2qmmm.py` into argument parsing (`cli.py`) and job execution (`app.py`) part
    * Removed `assess_job` function that did nothing
    * Formatted to a more conventional standard (as an example to discuss possible future improvements)
    * Used `pathlib` instead of `os`
  * Removed unnecessary shebangs (#!/usr/bin/env python)
  * Added directories for package structure `docs`, `examples`, `tests`, `lib`
  * Excluded `LibFiles` directory (What are these for? Should these be moved to examples or tests?)
  * Moved `json_files` to `lib/json`
  * Moved `correction_database` to `lib/correction_database`
  * Renamed modules and packages to be lower-case to meet the common convention
  * Excluded `Documentation` because there is nothing in there (and we should probably use the web page docs instead; where are these at the moment?)
  * Moved `Output` under `logging`
  * Changed super confusing imports like `sum_pcf_tm as make_p_sum`
