Developer notes
===============

Summary of changes in regards to the migration to v2
----------------------------------------------------

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
