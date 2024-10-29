import argparse
import pathlib
from datetime import datetime

import Generators.System as System
import Generators.GeneratorTopologies as Top
import Generators.GeneratorPCF as PCF
import Jobs.Singlepoint as SP
import Jobs.Optimization as OPT
from gmx2qmmm.types import StrPath
from Logging.Logger import Logger
from HandlingInput.HandlerInput import FileReader


#   TODO: Read And Include GROMACS / VMD Style Index Files; Extractable From Different Programs(?)
#   TODO: Include Prio For Library / Parameter Files

class App():
    """Main job orchestrator

    Args:
        parameters: Path to a file with input parameter definition
        logfile: Path to a logfile. If `None`, does not log any information.
        work_dir: Job working and output directory
    """

    def __init__(self, parameters: StrPath, logfile: StrPath | None = None, work_dir: StrPath = ".") -> None:

        # TODO: Work with parameters alternatively provided in a dictionary

        self.work_dir = pathlib.Path(work_dir).resolve()

        self.logfile = None
        if logfile is not None:
            self.logfile = pathlib.Path(logfile).resolve()

        log_entry ='Initialisation of gmx2qmmm completed'
        if (self.logfile is not None) and self.logfile.is_file():
            backup = self.backup_path(self.logfile)
            log_entry += f'; Backed up existing logfile to {backup}.'

        Logger.log_append(
            path_file_logging=self.logfile,
            str_info_to_log=log_entry
        )

        self.parameters = FileReader.read_file(
            path_file_logging=self.logfile,
            str_filename_input=parameters
        )

        log_entry = 'Parameters:\n' + '\n'.join([f'\t{k+ ":":<25}{v}' for k, v in self.parameters.items()])
        Logger.log_append(
            path_file_logging=self.logfile,
            str_info_to_log=log_entry
        )

        self.initalize_system()
        self.generate_topology()

        # Add list of elements to system instance
        self.system.list_atom_elements = self.system.get_atoms(self.topology.qmmm_topology)
        # NOTE: This operation (setting an instance attribute with the return value of a method on the very same instance) looks weird.
        #       Maybe that's a hint that the relationship between the System and Topology class is not ideal (JJ)

        self.generate_PCF()

    @classmethod
    def from_cli_args(cls, args: argparse.Namespace) -> 'App':
        """Initialise with CLI arguments parsed by argparse

        Args:
            args: An argparse namespace object obtained from a suitable
                parser
        """
        return cls(parameters=args.parameters, logfile=args.logfile)

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

    def initalize_system(self) -> None:
        """Setup system using input parameters"""
        self.system = System.SystemInfo(self.parameters)

    def generate_topology(self) -> None:
        """Setup topology using input parameters and system"""
        self.topology = Top.GenerateTopology(self.parameters, self.system, self.work_dir)

    def generate_PCF(self) -> None:
        """Setup point charge field using input parameters, system, and topology"""
        # NOTE (AJ): I forgot what I used the Job keyword for, I think I will only need it later, I will get back to that
        # NOTE (JJ): I find the note above very confusing. Can you please clarify or clean that up?
        self.pointchargefield = PCF.GeneratePCF(self.parameters, self.system, self.topology, self.work_dir)

    def run(self) -> None:
        """Run requested job"""

        job_func_map = {
            'singlepoint': SP.Singlepoint,
            'opt': OPT.Optimisation,
        }
        # TODO: Support NMA, scan, and opt_root_following (again)

        try:
            job_func_map[self.parameters['jobtype']](
                self.parameters, self.system, self.topology, self.pointchargefield, self.work_dir
                )
        except KeyError as exc:
            raise ValueError(
                f'Job type {self.parameters['jobtype']!r} not understood. '
                f'Must be one of {list(job_func_map.keys())}'
            ) from exc
