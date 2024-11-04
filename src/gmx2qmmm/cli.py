"""CLI interface"""

import argparse
import textwrap
from gmx2qmmm.app import App


def main():
    """Main entry point"""

    parser = argparse.ArgumentParser(
        description="gmx2qmmm, a Python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\n
            gmx2qmmm (GNU licenced) by Jan Philipp Goetze.\n
            gmx2qmmm is a free QM/MM interface for Gromacs. Compatible with most versions, starting from 5.X.\n
            gmx2qmmm uses an additive QM/MM scheme with a numerical charge shift procedure and correction of cut bonds via model potentials.\n
            """
        )
    )

    parser.add_argument(
        "-l",
        "--logfile",
        help="Logfile (.log)",
        type=str,
        default="logfile.txt"
    )

    parser.add_argument(
        "-p",
        "--parameters",
        help="Parameter file (.txt)",
        type=str,
        default="params.txt"
    )

    args = parser.parse_args()

    app = App.from_cli_args(args)
    app.run()
