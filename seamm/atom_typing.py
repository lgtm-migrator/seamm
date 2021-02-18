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
        self.atom_typer = None

    def create(self,
            typing_engine=None,
            forcefield=None,
            parameter_set=None
            ):


        self.typing_engine = typing_engine
        self.forcefield = forcefield 
        self.parameter_set = parameter_set 

        # Setup the plugin handling
        if self.typing_engine is None:
            raise NameError("Please provide an atom typing engine")

        if self.forcefield is None:
            raise NameError("Please provide a forcefield name")

        if self.parameter_set is None:
            raise NameError(f"Please provide a parameter set within the {self.forcefield} forcefield")

        atom_typer_plugin = self.plugin_manager.get(self.typing_engine)
        self.atom_typer = atom_typer_plugin.create_atom_typer(forcefield=self.forcefield, 
                parameter_set=self.parameter_set)
        return self.atom_typer
