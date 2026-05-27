import argparse
import pathlib
import shutil
from datetime import datetime
from typing import Optional

from loguru import logger

from gmx2qmmm.generators.system import SystemInfo
from gmx2qmmm.generators import topology
from gmx2qmmm.generators import pcf
from gmx2qmmm.jobs.singlepoint import Singlepoint
from gmx2qmmm.jobs.optimisation import Optimisation
from gmx2qmmm.types import StrPath
from gmx2qmmm.input_handler import FileReader


#   TODO: Read And Include GROMACS / VMD Style Index Files; Extractable From Different Programs(?)
#   TODO: Include Prio For Library / Parameter Files


class App:
    """Main job orchestrator

    Args:
        parameters: Path to a file with input parameter definition
        logfile: Path to a logfile. If `None`, use `"logfile.log"` in the `working_directory`.
        work_dir: Job working and output directory
    """

    def __init__(
        self,
        parameters: StrPath,
        logfile: Optional[StrPath] = None,
        loglevel: str = "DEBUG",
        work_dir: StrPath = ".",
    ) -> None:
        # TODO: Work with parameters alternatively provided in a dictionary
        self.base_dir = pathlib.Path(__file__).parents[2]
        self.work_dir = pathlib.Path(work_dir).resolve()

        self.loglevel = loglevel
        if logfile is None:
            logfile = self.work_dir / "logfile.log"
        self.logfile = logfile

        logger.info("Initialising gmx2qmmm")

        self.parameters = FileReader.read_file(
            path_file_logging=self.logfile, str_filename_input=parameters
        )
        self.check_gmx_paths()

        log_entry = "Parameters:\n" + "\n".join(
            [f"\t{k + ':':<25}{v}" for k, v in self.parameters.items()]
        )
        logger.info(log_entry)

        self.initialize_system()
        self.generate_topology()

        # Add list of elements to system instance
        self.system.list_atom_elements = self.system.get_atoms(
            self.topology.qmmm_topology
        )
        # NOTE: This operation (setting an instance attribute with the return value of a method on the very same instance) looks weird.
        #       Maybe that's a hint that the relationship between the System and Topology class is not ideal (JJ)

        self.generate_PCF()

    @property
    def logfile(self) -> pathlib.Path:
        return self._logfile

    @logfile.setter
    def logfile(self, value: StrPath) -> None:
        self._logfile = pathlib.Path(value).resolve()
        if self._logfile.is_file():
            backup = self.backup_path(self.logfile)
        else:
            backup = None

        logger.remove()
        logger.add(self._logfile, level=self.loglevel)

        if backup is not None:
            logger.info(f"Backed up existing logfile to {backup}")

    @classmethod
    def from_cli_args(cls, args: argparse.Namespace) -> "App":
        """Initialise with CLI arguments parsed by argparse

        Args:
            args: An argparse namespace object obtained from a suitable
                parser
        """
        return cls(parameters=args.parameters, logfile=args.logfile, work_dir=args.working_directory)

    @staticmethod
    def backup_path(path: pathlib.Path) -> pathlib.Path:
        """Backup a file or directory

        Includes the current timestamp in the backup path.

        Args:
            path: Path to backup

        Returns:
            Backup path
        """
        now = datetime.now()
        timestamp = f"{now.date()}_{now.hour}-{now.minute}-{now.second}"
        return path.rename(path.with_name(f"{path.name}_{timestamp}"))

    def initialize_system(self) -> None:
        """Setup system using input parameters"""
        self.system = SystemInfo(self.parameters, self.base_dir, self.work_dir)
        logger.success("Initialised system")

    def generate_topology(self) -> None:
        """Setup topology using input parameters and system"""
        self.topology = topology.GenerateTopology(
            self.parameters, self.system, self.work_dir
        )
        logger.success("Generated topology")

    def generate_PCF(self) -> None:
        """Setup point charge field using input parameters, system, and topology"""
        self.pointchargefield = pcf.GeneratePCF(
            self.parameters, self.system, self.topology, self.work_dir
        )
        logger.success("Generated PCF")

    def run(self) -> None:
        """Run requested job"""

        job_func_map = {
            "singlepoint": Singlepoint,
            "opt": Optimisation,
        }
        # TODO: Support NMA, scan, and opt_root_following (again)

        try:
            job_func = job_func_map[self.parameters["jobtype"]]
        except KeyError as exc:
            raise ValueError(
                f"Job type {self.parameters['jobtype']!r} not understood. "
                f"Must be one of {list(job_func_map.keys())}"
            ) from exc

        job_func(
            self.parameters,
            self.system,
            self.topology,
            self.pointchargefield,
            self.work_dir,
            self.base_dir,
        )

    def check_gmx_paths(self) -> None:
        """Set default paths for GROMACS executables and files if not provided in parameters"""

        self.parameters.setdefault("gmxcmd", "gmx")
        gmxcmd = self.parameters["gmxcmd"]
        gmx_path = shutil.which(gmxcmd)

        if "gmxpath" not in self.parameters or "gmxtop_path" not in self.parameters:

            if gmx_path is None:
                raise ValueError(
                    f"Could not find {gmxcmd!r} in PATH.\n"
                    "Please provide the path to the GROMACS executables"
                    "and/or bin/top paths in the parameters."
                    )

            if "gmxpath" not in self.parameters:
                self.parameters["gmxpath"] = str(pathlib.Path(gmx_path).parent)

            if "gmxtop_path" not in self.parameters:
                self.parameters["gmxtop_path"] = str(pathlib.Path(gmx_path).parent.parent / "share" / "gromacs" / "top")
