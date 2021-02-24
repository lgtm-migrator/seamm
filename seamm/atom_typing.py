# -*- coding: utf-8 -*-

import json
from datetime import datetime
import logging
import seamm
import seamm_util  # MUST come after seamm
import os
import os.path
import stat

logger = logging.getLogger(__name__)

class AtomTyperFactory:

    def __init__(
        self,
        namespace='org.molssi.seamm.atom_typers',
    ):
        self.plugin_manager = seamm.PluginManager(namespace)

    def create(self,
            atomtyping_engine=None,
            forcefield=None,
            parameter_set=None
            ):


        self.atomtyping_engine = atomtyping_engine
        self.forcefield = forcefield 
        self.parameter_set = parameter_set 

        # Setup the plugin handling
        if self.atomtyping_engine is None:
            raise NameError("Please provide an atom typing engine")

        if self.forcefield is None:
            raise NameError("Please provide a forcefield name")

        if self.parameter_set is None:
            raise NameError(f"Please provide a parameter set within the {self.forcefield} forcefield")

        atomtyper_plugin = self.plugin_manager.get(self.atomtyping_engine)
        self.atomtyper = atomtyper_plugin.create_atomtyper(forcefield=self.forcefield, 
                parameter_set=self.parameter_set)
        return self.atomtyper
