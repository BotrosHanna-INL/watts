from datetime import datetime
import os
from pathlib import Path
import shutil
import subprocess
import time

import h5py
import numpy as np
import pandas as pd

from .fileutils import PathLike
from .parameters import Parameters
from .plugin import TemplatePlugin
from .results import Results


class ResultsSAM(Results):
    def __init__(self, params: Parameters, time, inputs, outputs):
        super().__init__('SAM', params, time, inputs, outputs)

    @property
    def stdout(self):
        return (self.base_path / "SAM_log.txt").read_text()

    def save(self, filename: PathLike):
        """Save results to an HDF5 file

        Parameters
        ----------
        filename
            File to save results to
        """
        with h5py.File(filename, 'w') as h5file:
            super()._save(h5file)

    @classmethod
    def _from_hdf5(cls, obj: h5py.Group):
        """Load results from an HDF5 file

        Parameters
        ----------
        obj
            HDF5 group to load results from
        """
        time, parameters, inputs, outputs = Results._load(obj)
        return cls(parameters, time, inputs, outputs)


class PluginSAM(TemplatePlugin):
    """Plugin for running SAM

    Parameters
    ----------
    template_file
        Templated SAM input
    """
    def  __init__(self, template_file: str):
        super().__init__(template_file)
        self._sam_exec = 'sam-opt'

    @property
    def sam_exec(self):
        return self._sam_exec

    @sam_exec.setter
    def sam_exec(self, exe: PathLike):
        if shutil.which(exe) is None:
            raise RuntimeError(f"SAM executable '{exe}' is missing.")
        self._sam_exec = Path(exe)

    def options(self, sam_exec):
        """Input SAM user-specified options

        Parameters
        ----------
        SAM_exec
            Path to SAM executable
        """
        self.sam_exec = sam_exec
        self.sam_inp_name = "SAM.i"

    def prerun(self, model: Parameters):
        """Generate the SAM input based on the template

        Parameters
        ----------
        model
            Model used when rendering template

        """
        self._run_time = time.time_ns()
        # Render the template
        print("Pre-run for SAM Plugin")
        super().prerun(model, filename=self.sam_inp_name)

    def run(self):
        """Run SAM"""
        print("Run for SAM Plugin")

        log_file = Path("SAM_log.txt")

        # Run SAM and store  error message to SAM log file
        with log_file.open("w") as outfile:
            subprocess.run(
                [self.sam_exec, "-i", self.sam_inp_name],
                stdout=outfile,
                stderr=subprocess.STDOUT
            )

    def postrun(self, model: Parameters) -> ResultsSAM:
        """Read SAM results and store in model

        Parameters
        ----------
        model
            Model to store SAM results in
        """
        print("post-run for SAM Plugin")
        self._save_SAM_csv(model)

        time = datetime.fromtimestamp(self._run_time * 1e-9)
        inputs = ['SAM.i']
        outputs = [p for p in Path.cwd().iterdir() if p.name not in inputs]
        return ResultsSAM(model, time, inputs, outputs)

    def _save_SAM_csv(self, model):
        """Read all SAM '.csv' files and store in model

        Parameters
        ----------
        model
            Model to store SAM results in
        """
        csv_file_name = self.sam_inp_name[:-2] + "_csv.csv"
        # Save SAM's main output '.csv' files
        if os.path.isfile(csv_file_name):
            csv_file_df = pd.read_csv(csv_file_name)
            for column_name in csv_file_df.columns:
                model.set(column_name, np.array(csv_file_df[column_name]), user='plugin_sam')

        # Read SAM's vector postprocesssor '.csv' files and save the parameters as individual array
        exist_name = []
        for file in os.listdir():
            if file.startswith(self.sam_inp_name[:-2] + "_csv_") and not file.endswith("_0000.csv"):
                vector_csv_df = pd.read_csv(file)
                csv_param = list(set(list(vector_csv_df.columns)) - set(set(["id", "x", "y", "z"])))
                model.set(file[:-4], np.array(vector_csv_df[csv_param[0]]).astype(np.float64), user='plugin_sam')

                for name in ["id", "x", "y", "z"]:
                    new_name = file[:-8] + name
                    if new_name not in exist_name:
                        model.set(new_name, np.array(vector_csv_df[name]).astype(np.float64), user='plugin_sam')
                        exist_name.append(file[:-8] + name)
