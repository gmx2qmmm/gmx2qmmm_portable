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
  * Removed unnecessary shebangs (#!/usr/bin/env python)
  * Started to use `pathlib` instead of `os`