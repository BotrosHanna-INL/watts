# SPDX-FileCopyrightText: 2022 UChicago Argonne, LLC
# SPDX-License-Identifier: MIT

import os
import glob
import sys
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pathlib import Path
import shutil
import time
from typing import List, Optional

import json
import csv
import numpy as np
import pandas as pd
import subprocess

from .fileutils import PathLike, run as run_proc, tee_stdout, tee_stderr
from .parameters import Parameters
from .plugin import TemplatePlugin
from .results import Results


class ResultsDakota(Results):
    """Dakota simulation results

    Parameters
    ----------
    params
        Parameters used to generate inputs
    name
        Name of workflow producing results
    time
        Time at which workflow was run
    inputs
        List of input files
    outputs
        List of output files

    Attributes
    ----------
    stdout
        Standard output from Dakota run
    output_data
        Dictionary with data from .csv files
    """
    def __init__(self, params: Parameters, name: str, time: datetime,
                 inputs: List[PathLike], outputs: List[PathLike]):
        super().__init__('Dakota', params, name, time, inputs, outputs)
        self.output_data = self._get_Dakota_csv()

    @property
    def stdout(self) -> str:
        return (self.base_path / "Dakota_log.txt").read_text()

    def _get_Dakota_csv(self) -> dict:
        """Read all Dakota '.dat' files and return results in a dictionary

        Returns
        -------
        Results from Dakota .dat files

        """

        # Save Dakota's main output '.dat' files
        output_data = {}
        
        if (glob.glob('dakota_opt.dat')):
            with open('dakota_opt.dat') as f:
                reader=csv.reader(f)
                rows=[row for idx, row in enumerate(reader) if idx == 0]
      
            col_names = str(rows[0][0]).split()
            df = pd.read_csv("dakota_opt.dat", sep="\s+", skiprows=1, names=col_names)

            for name in col_names:
                output_data[name] = np.array(df[name])

        # Save Dakota's final output '.dat' files
        if (glob.glob('finaldata1.dat')):
            with open('finaldata1.dat') as fd:
                reader=csv.reader(fd)
                rows=[row for idx, row in enumerate(reader) if idx == 0]

            new_rows = str(rows[0][0]).split()
            output_data['finaldata1'] = np.array([float(i) for i in new_rows])

        return output_data


class PluginDakota(TemplatePlugin):
    """Plugin for running Dakota

    Parameters
    ----------
    template_file
        Templated Dakota input
    show_stdout
        Whether to display output from stdout when Dakota is run
    show_stderr
        Whether to display output from stderr when Dakota is run
    extra_inputs
        List of extra (non-templated) input files that are needed
    extra_template_inputs
        Extra templated input files

    Attributes
    ----------
    dakota_exec
        Path to Dakota executable

    """

    def __init__(self, template_file: str, show_stdout: bool = False,
                 show_stderr: bool = False, 
                 extra_inputs: Optional[List[str]] = None,
                 extra_template_inputs: Optional[List[PathLike]] = None):
        super().__init__(template_file, extra_inputs, extra_template_inputs)

        dakota_dir = Path(os.environ.get("DAKOTA_DIR", ""))
        self._dakota_exec = dakota_dir / f"dakota.sh"

        self.dakota_inp_name = "dakota_watts_opt.in"
        self.show_stdout = show_stdout
        self.show_stderr = show_stderr

        self._initial_dir = os.getcwd()

    @property
    def dakota_exec(self) -> Path:
        return self._dakota_exec

    @dakota_exec.setter
    def dakota_exec(self, exe: PathLike):
        if shutil.which(exe) is None:
            raise RuntimeError(f"Dakota executable '{exe}' is missing.")
        self._dakota_exec = Path(exe)

    def options(self, dakota_exec):
        """Input Dakota user-specified options

        Parameters
        ----------
        Dakota_exec
            Path to DAKOTA executable
        """
        self.dakota_exec = dakota_exec

    def prerun(self, params: Parameters):
        """Generate the Dakota input based on the template

        Parameters
        ----------
        params
            Parameters used when rendering template
        """
        # Render the template
        # Make a copy of params and convert units if necessary
        # The original params remains unchanged

        params_copy = params.convert_units()

        print("Pre-run for Dakota Plugin")
        self._run_time = time.time_ns()
        super().prerun(params_copy, filename=self.dakota_inp_name)

        # Copy all files to the temporary directory.
        # Avoid duplicate files.

        files = os.listdir(self._initial_dir)
        _current_dir = os.getcwd()
        for fname in files:
            if not os.path.isdir(os.path.join(self._initial_dir, fname)):
                if not os.path.exists(os.path.join(_current_dir, fname)):
                    shutil.copy2(os.path.join(self._initial_dir, fname), _current_dir)
                else:
                    print(os.path.join(self._initial_dir, fname), _current_dir + " already exists")   

    def run(self):
        """Run Dakota"""
        print("Run for Dakota Plugin")

        log_file = Path("Dakota.txt")

        with log_file.open("w") as outfile:
            func_stdout = tee_stdout if self.show_stdout else redirect_stdout
            func_stderr = tee_stderr if self.show_stderr else redirect_stderr
            with func_stdout(outfile), func_stderr(outfile):
                run_proc([self.dakota_exec, "-i", self.dakota_inp_name])

    def postrun(self, params: Parameters, name: str) -> ResultsDakota:
        """Read Dakota results and create results object

        Parameters
        ----------
        params
            Parameters used to create Dakota model
        name
            Name of the workflow

        Returns
        -------
        Dakota results object
        """
        print("Post-run for Dakota Plugin")

        time, inputs, outputs = self._get_result_input(self.dakota_inp_name)
        return ResultsDakota(params, name, time, inputs, outputs)

class PluginDakotaDriver():
    """Plugin for running the Dakota driver

    Parameters
    ----------
    Attributes
    ----------
    """

    def __init__(self, coupled_code_exec: str):

        self.coupled_code_exec = coupled_code_exec

        self.parse_dakota_input()
        self.run_coupled_code()
        self.return_dakota_input()

    def parse_dakota_input(self):
        """Parse Dakota input

        Parameters
        ----------

        Returns
        -------
        """
        from interfacing import interfacing as di # Dakota's interface module        

        params, self.results = di.read_parameters_file()

        # Dump params to external params.json file for future use by the template engine
        params_for_template_engine_file_path = "params.json"
        with open(params_for_template_engine_file_path, 'w') as outfile:
            f = json.dump(params._variables,  outfile, default=lambda o: o.__dict__)
        
        f = open("params.json")
        data = json.load(f)
        f.close()

        f = open("input.txt", "w+")
        with open('input.tmpl') as file:
            while (line := file.readline().rstrip()):

                i_0 = line.index('=')
                i_1 = line.index('<') + 1
                i_2 = line.index('>') 

                f.write(line[0:i_0] + " = " + str(data[line[i_1:i_2]]) + "\n")
        f.close()

    def run_coupled_code(self):
        if os.path.exists(self.coupled_code_exec):
            P = subprocess.check_output(["python", self.coupled_code_exec])
        else:
            raise RuntimeError("Coupled-code script missing.")

        res_output = []
        with open('opt_res.out') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) > 0:
                    new_str = str(row).split()
                    for j in new_str:
                        try:
                            res_output.append(float(j))
                        except:
                            pass

        self.retval = dict([])
        self.retval['fns'] = res_output

    def return_dakota_input(self):

        # Insert extracted values into results
        # Results iterator provides an index, response name, and response
        try:
          for i, n, r in self.results:
            if r.asv.function:
                try:
                    r.function = self.retval['fns'][i]
                except:
                    pass
        # Catch Dakota 6.9 exception where results interface has changed
        # ValueError: too many values to unpack
        except ValueError:
          i = 0
          for n, r in self.results.items():
              r.function = self.retval['fns'][i]
              i+=1
        self.results.write()

        # Dump to external results.json file
        with open('results.json', 'w') as outfile:
            rst = json.dump(self.results,  outfile, default=lambda o: o.__dict__)

