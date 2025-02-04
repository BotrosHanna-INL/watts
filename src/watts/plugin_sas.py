# SPDX-FileCopyrightText: 2022-2023 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

import glob
import os
from pathlib import Path
import shutil
import subprocess
from typing import List, Optional

import numpy as np
import pandas as pd

from .fileutils import PathLike
from .parameters import Parameters
from .plugin import PluginGeneric, _find_executable
from .results import Results, ExecInfo


class ResultsSAS(Results):
    """SAS simulation results

    Parameters
    ----------
    params
        Parameters used to generate inputs
    exec_info
        Execution information (job ID, plugin name, time, etc.)
    inputs
        List of input files
    outputs
        List of output files

    Attributes
    ----------
    stdout
        Standard output from SAS run
    csv_data
        Dictionary with data from .csv files
    """
    def __init__(self, params: Parameters, exec_info: ExecInfo,
                 inputs: List[PathLike], outputs: List[PathLike]):
        super().__init__(params, exec_info, inputs, outputs)
        self.csv_data = self._get_sas_csv_data()

    def _get_sas_csv_data(self) -> dict:
        """Read all sas '.csv' files and return results in a dictionary

        Returns
        -------
        Results from sas .csv files

        """
        csv_files = glob.glob("*.csv") # List all csv files

        csv_data = {}
        for file in csv_files:
            if os.path.getsize(file) > 0: # Check if file is empty
                csv_data_sub = {}
                csv_file_df = pd.read_csv(file)

                for column_name in csv_file_df.columns:
                    csv_data_sub[column_name] =  np.array(csv_file_df[column_name])

                csv_data[file[:-4]] = csv_data_sub # Save sub dictionary as csv file name
                                                   # Removed .csv extension from file name

        return csv_data


class PluginSAS(PluginGeneric):
    """Plugin for running SAS

    Parameters
    ----------
    template_file
        Templated SAS input
    executable
        Path to SAS executable
    extra_inputs
        List of extra (non-templated) input files that are needed
    extra_template_inputs
        Extra templated input files
    show_stdout
        Whether to display output from stdout when SAS is run
    show_stderr
        Whether to display output from stderr when SAS is run

    Attributes
    ----------
    executable
        Path to SAS executable
    execute_command
        List of command-line arguments used to call the executable
    conv_channel
        Path to CHANNELtoCSV utility executable
    conv_mfuel
        Path to MFUELPeaktoCSV utility executable
    conv_primar4
        Path to PRIMAR4toCSV utility executable

    """

    def  __init__(
        self,
        template_file: str,
        executable: PathLike = 'sas.x',
        extra_inputs: Optional[List[str]] = None,
        extra_template_inputs: Optional[List[PathLike]] = None,
        show_stdout: bool = False,
        show_stderr: bool = False
    ):
        executable = _find_executable(executable, 'SAS_DIR')
        execute_command = ['{self.executable}', '-i', '{self.input_name}',
                           '-o', 'out.txt']
        super().__init__(executable, execute_command, template_file, extra_inputs,
                         extra_template_inputs, "SAS", show_stdout, show_stderr)
        self.input_name = "SAS.inp"

        # Set other executables based on the main SAS executable
        suffix = executable.suffix
        self._conv_channel = executable.with_name(f"CHANNELtoCSV{suffix}")
        self._conv_mfuel = executable.with_name(f"MFUELtoCSV{suffix}")
        self._conv_primar4 = executable.with_name(f"PRIMAR4toCSV{suffix}")

    @property
    def conv_channel(self) -> Path:
        return self._conv_channel

    @property
    def conv_primar4(self) -> Path:
        return self._conv_primar4

    @property
    def conv_mfuel(self) -> Path:
        return self._conv_mfuel

    @conv_channel.setter
    def conv_channel(self, exe: PathLike):
        if shutil.which(exe) is None:
            raise RuntimeError(f"CHANNELtoCSV utility executable '{exe}' is missing.")
        self._conv_channel = Path(exe)

    @conv_primar4.setter
    def conv_primar4(self, exe: PathLike):
        if shutil.which(exe) is None:
            raise RuntimeError(f"PRIMAR4toCSV utility executable '{exe}' is missing.")
        self._conv_primar4 = Path(exe)

    @conv_mfuel.setter
    def conv_mfuel(self, exe: PathLike):
        if shutil.which(exe) is None:
            raise RuntimeError(f"MFUELtoCSV utility executable '{exe}' is missing.")
        self._conv_mfuel = Path(exe)

    @property
    def execute_command(self):
        return [str(self.executable), "-i", self.input_name, "-o", "out.txt"]

    def postrun(self, params: Parameters, exec_info: ExecInfo) -> ResultsSAS:
        """Read SAS results and create results object

        Parameters
        ----------
        params
            Parameters used to create SAS model
        name
            Name of the workflow

        Returnss
        -------
        SAS results object
        """

        # Convert CHANNEl.dat and PRIMAR4.dat and MFUELss.dat to csv files
        # using SAS utilities. Check if files exist because
        # they may not be outputted per user's choice.

        if Path("CHANNEL.dat").is_file():
            with open("CHANNEL.dat", "r") as file_in, open("CHANNEL.csv", "w") as file_out:
                subprocess.run(str(self.conv_channel), stdin=file_in, stdout=file_out)

        if Path("PRIMAR4.dat").is_file():
            with open("PRIMAR4.dat", "r") as file_in, open("PRIMAR4.csv", "w") as file_out:
                subprocess.run(str(self.conv_primar4), stdin=file_in, stdout=file_out)

        if Path("MFUELss_C000001.dat").is_file(): 
            with open("MFUELss_C000001.dat", "r") as file_in: 
                subprocess.run(str(self.conv_mfuel), stdin=file_in) 
   
        return super().postrun(params, exec_info)
