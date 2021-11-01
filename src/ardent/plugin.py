from abc import ABC, abstractmethod
import os, sys, shutil
import csv


class Plugin(ABC):
    def __init__(self, template_file):
        ...

    @abstractmethod
    def prerun(self, model):
        # Fill in the template to create real inputs

        # Run arbitrary user scripts
        ...
        
    @abstractmethod
    def prerun(self, model):
        # Fill in the template to create real inputs

        # Run arbitrary user scripts
        ...

    @abstractmethod
    def run(self):
        ...

    @abstractmethod
    def postrun(self, model):
        ...


class OpenmcPlugin(Plugin):
    def __init__(self, model_builder):
        self.model_builder = model_builder

    def prerun(self, model):
        self.model_builder(model)


class ExamplePlugin(Plugin):
    def  __init__(self, model_builder):
        self.model_builder = model_builder

    def workflow(self, model):
        prerun_crash = self.prerun(model)
        run_crash = self.run()
        postrun_crash = self.postrun()

    def prerun(self, model):
        # Render the template
        print("Pre-run for Example Plugin")
        self.model_builder(model)

    def run(self):
        print("Run for Example Plugin")
        

    def postrun(self):
        print("post-run for Example Plugin")

class PluginSAM(Plugin):
    def  __init__(self, model_builder):
        self.model_builder = model_builder
        self.sam_inp_name = "SAM.i"
        self.sam_tmp_folder = "tmp_SAM"
        self.SAM_exec = "../sam-opt-mpi"

    def workflow(self, model):
        prerun_crash = self.prerun(model)
        run_crash = self.run()
        postrun_crash = self.postrun(model)

    def prerun(self, model):
        # Render the template
        print("Pre-run for SAM Plugin")
        self.model_builder(model)

    def run(self):
        print("Run for SAM Plugin")

        if os.path.exists(self.sam_tmp_folder):
            shutil.rmtree(self.sam_tmp_folder)
        os.mkdir(self.sam_tmp_folder)

        shutil.copy("sam_template.rendered", self.sam_tmp_folder+"/"+self.sam_inp_name)
        os.chdir(self.sam_tmp_folder)
        
        if os.path.isfile(self.SAM_exec) is False:
            raise RuntimeError("SAM executable missing")
        os.system(self.SAM_exec+" -i "+self.sam_inp_name+" > "+self.sam_inp_name+".out")

    def postrun(self, model):
        print("post-run for SAM Plugin")
        # TODO: find all '.cvs' files
        csv_file_name = "SAM_csv.csv"

        if os.path.isfile(csv_file_name):
            cvs_file = open(csv_file_name,'r')
            read_file = csv.reader(cvs_file)
            out_press = None
            for row in read_file:
                # TODO: save all results from CSV files - this example only save 1 result
                if row[0] == "1":
                    model['max_Pcoolant'] = float(row[1])

        os.chdir("../")
        shutil.rmtree(self.sam_tmp_folder)

