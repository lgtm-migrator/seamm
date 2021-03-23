import logging
logger = logging.getLogger(__name__)

class EnergyExpression:


    def __init__(self, system=None, atomtyping_engine=None, configuration=None, style=''):
        """Create the energy expression for the given structure

        Parameters
        ----------
        system : _System
            The system being used
        configuration : int = None
            Which configuration. Defaults to the current_configuration.
        style : str = ''
            The style of energy expression. Currently only 'LAMMPS' is
            supported.

        Returns
        -------
        self.eex : {str: []}
            The energy expression as a dictionary of terms
        """
        logger.debug('Creating the self.eex')

        self.eex = {}

        self.system = system

        self.atomtyping_engine = atomtyping_engine 
        self.current_forcefield = self.atomtyping_engine.name

        if configuration is None:
            configuration = system.current_configuration

        # We will need the elements for fix shake, 1-based.
        # self.eex['elements'] = ['']
        # self.eex['elements'].extend(self.system['atom'].symbols(configuration))

        # The periodicity & cell parameters
        # self.eex['periodicity'] = self.system.periodicity
        self.eex['periodicity'] = self.system.system.configuration.periodicity

        if self.eex['periodicity'] == 3:
            # self.eex['cell'] = self.system.system['cell'].cell().parameters
            self.eex['cell'] = self.system.system.configuration.cell.parameters

        self.setup_topology(system, configuration, style)

        self.eex_atoms(self.eex, system, configuration)

        logger.debug(f'    forcefield terms: {self.atomtyping_engine.forcefield.ff["terms"]}')

        for term in self.atomtyping_engine.forcefield.ff['terms']:
            function_name = 'eex_' + term.replace('-', '_')
            function_name = function_name.replace(' ', '_')
            function_name = function_name.replace(',', '_')
            function = getattr(self, function_name, None)
            if function is None:
                print('Function {} does not exist yet'.format(function_name))
            else:
                function(self.eex, system, configuration)

        # Now run through the sections for the functionals forms,
        # processing each
        for fform in self.atomtyping_engine.forcefield.ff['functional_forms']:
            self.atomtyping_engine.forcefield._get_parameters(fform, self.atomtyping_engine.forcefield.version)

        if logger.isEnabledFor(logging.DEBUG):
            section = 'bond_increments'
            try:
                print(json.dumps(self.atomtyping_engine.forcefield.ff[section], indent=4))
            except:  # noqa: E722
                pprint.pprint(self.atomtyping_engine.forcefield.ff[section])

    def setup_topology(self, system, configuration=None, style=''):
        """Create the list of bonds, angle, torsion, etc. for the system

        This topology information is held in self.topology.

        Parameters
        ----------
        system : _System
            The system being used
        configuration : int = None
            Which configuration. Defaults to the current_configuration.
        style : str = ''
            The style of energy expression. Currently only 'LAMMPS' is
            supported.

        Returns
        -------
        None
        """
        self.topology = {}

        if configuration is None:
            configuration = system.current_configuration

        sys_atoms = system['atom']
        sys_bonds = system['bond']

        #n_atoms = sys_atoms.n_atoms(configuration)
        n_atoms = system.system.configuration.n_atoms
        self.topology['n_atoms'] = n_atoms

        # Need the transformation from atom ids to indices
        #atom_ids = sys_atoms.atom_ids(configuration)
        atom_ids = system.system.configuration.atoms.ids
        to_index = {j: i + 1 for i, j in enumerate(atom_ids)}

        # extend types with a blank so can use 1-based indexing
        types = self.topology['types'] = ['']
        key = f'atom_types_{self.atomtyping_engine.name}'
        #types.extend(sys_atoms.get_column(key, configuration=configuration))
        types.extend(sys_atoms.get_column(key))
        # bonds
        # bonds = self.topology['bonds'] = [
        #    (to_index[row['i']], to_index[row['j']])
        #    for row in sys_bonds.bonds(configuration=configuration)
        #]

        bonds = self.topology['bonds'] = [(x, y) for x, y in zip(sys_bonds.get_as_dict()['i'], sys_bonds.get_as_dict()['j'])]

        # atoms bonded to each atom i
        self.topology['bonds_from_atom'] = system.system.configuration.bonded_neighbors()
        bonds_from_atom = self.topology['bonds_from_atom']

        # angles
        angles = self.topology['angles'] = []

       # for j in range(1, n_atoms + 1):
       #     for i in bonds_from_atom[j]:
       #         for k in bonds_from_atom[j]:
       #             if i < k:
       #                 angles.append((i, j, k))

        for j in bonds_from_atom.keys():
            for i in bonds_from_atom[j]:
                for k in bonds_from_atom[j]:
                    if i < k:
                        angles.append((i, j, k))

        # torsions
        torsions = self.topology['torsions'] = []
        for j, k in bonds:
            for i in bonds_from_atom[j]:
                if i == k:
                    continue
                for l in bonds_from_atom[k]:  # noqa: E741
                    if l == j:  # noqa: E741
                        continue
                    torsions.append((i, j, k, l))

        # Out-of-planes
        oops = self.topology['oops'] = []
        for m in bonds_from_atom.keys():
            if len(bonds_from_atom[m]) == 3:
                i, j, k = bonds_from_atom[m]
                oops.append((i, m, j, k))

        if style == 'LAMMPS':
            for m in bonds_from_atom.keys():
                if len(bonds_from_atom[m]) == 4:
                    i, j, k, l = bonds_from_atom[m]  # noqa: E741
                    oops.append((i, m, j, k))
                    oops.append((i, m, j, l))
                    oops.append((i, m, k, l))
                    oops.append((j, m, k, l))

    def eex_atoms(self, eex, system, configuration=None):
        """List the atoms into the energy expression"""

        #atoms = system['atom']
        #coordinates = atoms.coordinates(
        #    configuration=configuration, fractionals=False
        #)
        atoms = system.systems[0].configuration.atoms
        coordinates = atoms.coordinates

        key = f'atomtypes_{self.current_forcefield}'
        types = atoms.get_column(key)

        result = eex['atoms'] = []
        atom_types = eex['atom types'] = []
        masses = eex['masses'] = []

        for itype, xyz in zip(types, coordinates):
            if itype in atom_types:
                index = atom_types.index(itype) + 1
            else:
                atom_types.append(itype)
                index = len(atom_types)
                masses.append((self.mass(itype), itype))
            x, y, z = xyz
            result.append((x, y, z, index))

        eex['n_atoms'] = len(result)
        eex['n_atom_types'] = len(atom_types)

    def eex_bond(self, eex, system, configuration=None):
        """Create the bond portion of the energy expression"""
        types = self.topology['types']
        bonds = self.topology['bonds']

        result = eex['bonds'] = []
        parameters = eex['bond parameters'] = []
        for i, j in bonds:
            parameters_type, real_types, form, parameter_values = \
                self.bond_parameters(types[i], types[j])
            new_value = (
                form, parameter_values, (types[i], types[j]), parameters_type,
                real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, index))
        eex['n_bonds'] = len(result)
        eex['n_bond_types'] = len(parameters)

    def eex_angle(self, eex, system, configuration=None):
        """Create the angle portion of the energy expression"""
        types = self.topology['types']
        angles = self.topology['angles']

        result = eex['angles'] = []
        parameters = eex['angle parameters'] = []
        for i, j, k in angles:
            parameters_type, real_types, form, parameter_values = \
                self.angle_parameters(types[i], types[j], types[k])
            new_value = (
                form, parameter_values, (types[i], types[j], types[k]),
                parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, index))
        eex['n_angles'] = len(result)
        eex['n_angle_types'] = len(parameters)

    def eex_torsion(self, eex, system, configuration=None):
        """Create the torsion portion of the energy expression"""
        types = self.topology['types']
        torsions = self.topology['torsions']

        result = eex['torsions'] = []
        parameters = eex['torsion parameters'] = []
        for i, j, k, l in torsions:
            parameters_type, real_types, form, parameter_values = \
                self.torsion_parameters(types[i], types[j], types[k], types[l])
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_torsions'] = len(result)
        eex['n_torsion_types'] = len(parameters)

    def eex_out_of_plane(self, eex, system, configuration=None):
        """Create the out-of-plane portion of the energy expression"""
        types = self.topology['types']
        oops = self.topology['oops']

        result = eex['oops'] = []
        parameters = eex['oop parameters'] = []
        for i, j, k, l in oops:
            parameters_type, real_types, form, parameter_values = \
                self.oop_parameters(types[i], types[j], types[k], types[l],
                                    zero=True)
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_oops'] = len(result)
        eex['n_oop_types'] = len(parameters)

    def eex_bond_bond(self, eex, system, configuration=None):
        """Create the bond-bond portion of the energy expression"""
        types = self.topology['types']
        angles = self.topology['angles']

        result = eex['bond-bond'] = []
        parameters = eex['bond-bond parameters'] = []
        for i, j, k in angles:
            parameters_type, real_types, form, parameter_values = \
                self.bond_bond_parameters(
                    types[i], types[j], types[k], zero=True)
            new_value = (
                form, parameter_values, (types[i], types[j], types[k]),
                parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, index))
        eex['n_bond-bond'] = len(result)
        eex['n_bond-bond_types'] = len(parameters)

    def eex_bond_angle(self, eex, system, configuration=None):
        """Create the bond-angle portion of the energy expression"""
        types = self.topology['types']
        angles = self.topology['angles']

        result = eex['bond-angle'] = []
        parameters = eex['bond-angle parameters'] = []
        for i, j, k in angles:
            parameters_type, real_types, form, parameter_values = \
                self.bond_angle_parameters(
                    types[i], types[j], types[k], zero=True)
            new_value = (
                form, parameter_values, (types[i], types[j], types[k]),
                parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, index))
        eex['n_bond-angle'] = len(result)
        eex['n_bond-angle_types'] = len(parameters)

    def eex_torsion_middle_bond(self, eex, system, configuration=None):
        """Create the middle_bond-torsion portion of the energy expression"""
        types = self.topology['types']
        torsions = self.topology['torsions']

        result = eex['middle_bond-torsion_3'] = []
        parameters = eex['middle_bond-torsion_3 parameters'] = []
        for i, j, k, l in torsions:
            parameters_type, real_types, form, parameter_values = \
                self.middle_bond_torsion_3_parameters(
                    types[i], types[j], types[k], types[l], zero=True)
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_middle_bond-torsion_3'] = len(result)
        eex['n_middle_bond-torsion_3_types'] = len(parameters)

    def eex_torsion_end_bond(self, eex, system, configuration=None):
        """Create the end_bond-torsion portion of the energy expression"""
        types = self.topology['types']
        torsions = self.topology['torsions']

        result = eex['end_bond-torsion_3'] = []
        parameters = eex['end_bond-torsion_3 parameters'] = []
        for i, j, k, l in torsions:
            parameters_type, real_types, form, parameter_values = \
                self.end_bond_torsion_3_parameters(
                    types[i], types[j], types[k], types[l], zero=True)
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_end_bond-torsion_3'] = len(result)
        eex['n_end_bond-torsion_3_types'] = len(parameters)

    def eex_torsion_angle(self, eex, system, configuration=None):
        """Create the angle-torsion portion of the energy expression"""
        types = self.topology['types']
        torsions = self.topology['torsions']

        result = eex['angle-torsion_3'] = []
        parameters = eex['angle-torsion_3 parameters'] = []
        for i, j, k, l in torsions:
            parameters_type, real_types, form, parameter_values = \
                self.angle_torsion_3_parameters(
                    types[i], types[j], types[k], types[l], zero=True)
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_angle-torsion_3'] = len(result)
        eex['n_angle-torsion_3_types'] = len(parameters)

    def eex_angle_torsion_angle(self, eex, system, configuration=None):
        """Create the angle-angle-torsion portion of the energy expression"""
        types = self.topology['types']
        torsions = self.topology['torsions']

        result = eex['angle-angle-torsion_1'] = []
        parameters = eex['angle-angle-torsion_1 parameters'] = []
        for i, j, k, l in torsions:
            parameters_type, real_types, form, parameter_values = \
                self.angle_angle_torsion_1_parameters(
                    types[i], types[j], types[k], types[l], zero=True)
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_angle-angle-torsion_1'] = len(result)
        eex['n_angle-angle-torsion_1_types'] = len(parameters)

    def eex_1_3_bond_bond(self, eex, system, configuration=None):
        """Create the 1,3 bond-bond portion of the energy expression"""
        types = self.topology['types']
        torsions = self.topology['torsions']

        result = eex['bond-bond_1_3'] = []
        parameters = eex['bond-bond_1_3 parameters'] = []
        for i, j, k, l in torsions:
            parameters_type, real_types, form, parameter_values = \
                self.bond_bond_1_3_parameters(
                    types[i], types[j], types[k], types[l], zero=True)
            new_value = (
                form, parameter_values,
                (types[i], types[j], types[k],
                 types[l]), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_bond-bond_1_3'] = len(result)
        eex['n_bond-bond_1_3_types'] = len(parameters)

    def eex_angle_angle(self, eex, system, configuration=None):
        """Create the angle-angle portion of the energy expression

        j is the vertex atom of the angles. For the angle-angle parameters
        the bond j-k is the common bond, i.e. the angles are i-j-k and j-k l
        """
        types = self.topology['types']
        oops = self.topology['oops']

        result = eex['angle-angle'] = []
        parameters = eex['angle-angle parameters'] = []
        for i, j, k, l in oops:
            parameters_type, real_types, form, parameter_values = \
                self.angle_angle_parameters(
                    types[i], types[j], types[k], types[l], zero=True)
            K1 = parameter_values['K']
            Theta10 = parameter_values['Theta10']
            Theta30 = parameter_values['Theta20']
            tmp = self.angle_angle_parameters(
                types[k], types[j], types[i], types[l], zero=True
            )[3]
            K2 = tmp['K']
            Theta20 = tmp['Theta20']
            tmp = self.angle_angle_parameters(
                types[i], types[j], types[l], types[k], zero=True
            )[3]
            K3 = tmp['K']
            new_value = (
                form, {
                    'K1': K1,
                    'K2': K2,
                    'K3': K3,
                    'Theta10': Theta10,
                    'Theta20': Theta20,
                    'Theta30': Theta30
                }, (types[i], types[j], types[k], types[l]), parameters_type,
                real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append((i, j, k, l, index))
        eex['n_angle-angle'] = len(result)
        eex['n_angle-angle_types'] = len(parameters)

    def mass(self, i):
        """Return the atomic mass for an atom type i"""
        if i in self.atomtyping_engine.forcefield.ff['atom_types']:
            return self.atomtyping_engine.forcefield.ff['atom_types'][i]['mass'] 

        raise RuntimeError('no atom type data for {}'.format(i))
