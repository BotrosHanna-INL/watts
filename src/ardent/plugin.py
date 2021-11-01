from abc import ABC, abstractmethod


class Plugin(ABC):
    """Class defining the Plugin interface"""

    def __init__(self, template_file):
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
    """Plugin for running OpenMC"""

    def __init__(self, model_builder):
        self.model_builder = model_builder

    def prerun(self, model):
        self.model_builder(model)


class ExamplePlugin(Plugin):
    def  __init__(self, model_builder):
        self.model_builder = model_builder

    def prerun(self, model):
        # Render the template
        self.model_builder(model)

    def run(self):
        ...

    def postrun(self):
        ...
