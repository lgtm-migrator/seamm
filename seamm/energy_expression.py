import logging
import pprint
import json
from copy import deepcopy
logger = logging.getLogger(__name__)

class EnergyExpression:


    def __init__(self, configuration=None, atomtyping_engine=None, style=''):
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

        self.atomtyping_engine = atomtyping_engine 

        self.current_forcefield = self.atomtyping_engine.name

        if configuration is None:
            raise TypeError("A configuration must be provided to the energy expression initalization method") 

        self.configuration = configuration

        self.eex['functional_forms'] = deepcopy(self.atomtyping_engine.forcefield.ff["terms"])

        self.eex['elements'] = configuration.atoms.symbols
        # We will need the elements for fix shake, 1-based.
        # self.eex['elements'] = ['']
        # self.eex['elements'].extend(self.system['atom'].symbols(configuration))

        # The periodicity & cell parameters
        # self.eex['periodicity'] = self.system.periodicity
        self.eex['periodicity'] = configuration.periodicity

        if self.eex['periodicity'] == 3:
            # self.eex['cell'] = self.system.system['cell'].cell().parameters
            self.eex['cell'] = configuration.cell.parameters

        self.setup_topology(configuration, style)

        self.eex_atoms(self.eex)

        logger.debug(f'    forcefield terms: {self.atomtyping_engine.forcefield.ff["terms"]}')

        for term in self.atomtyping_engine.forcefield.ff['terms']:
            function_name = 'eex_' + term.replace('-', '_')
            function_name = function_name.replace(' ', '_')
            function_name = function_name.replace(',', '_')
            function = getattr(self, function_name, None)
            if function is None:
                print('Function {} does not exist yet'.format(function_name))
            else:
                function(self.eex)

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

    def setup_topology(self, configuration=None, style=''):
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
            raise Exception 

        sys_atoms = configuration.atoms
        sys_bonds = configuration.bonds

        #n_atoms = sys_atoms.n_atoms(configuration)
        n_atoms = configuration.n_atoms
        self.topology['n_atoms'] = n_atoms

        # Need the transformation from atom ids to indices
        #atom_ids = sys_atoms.atom_ids(configuration)
        atom_ids = configuration.atoms.ids
        to_index = {j: i + 1 for i, j in enumerate(atom_ids)}

        # extend types with a blank so can use 1-based indexing
        key = f'atomtypes_{self.atomtyping_engine.forcefield.name}'
        sys_atoms_dict = sys_atoms.get_as_dict() 
        types = self.topology['types'] = {x: y for x, y in zip(sys_atoms_dict['id'], sys_atoms_dict[key])}
        #types.extend(sys_atoms.get_column(key, configuration=configuration))
        # bonds
        # bonds = self.topology['bonds'] = [
        #    (to_index[row['i']], to_index[row['j']])
        #    for row in sys_bonds.bonds(configuration=configuration)
        #]

        bonds = self.topology['bonds'] = [(x, y) for x, y in zip(sys_bonds.get_as_dict()['i'], sys_bonds.get_as_dict()['j'])]

        # atoms bonded to each atom i
        self.topology['bonds_from_atom'] = configuration.bonded_neighbors()
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

    def eex_atoms(self, eex):
        """List the atoms into the energy expression"""

        #atoms = system['atom']
        #coordinates = atoms.coordinates(
        #    configuration=configuration, fractionals=False
        #)
        atoms = self.configuration.atoms
        coordinates = atoms.coordinates

        types = self.topology['types']

        result = eex['atoms'] = []
        atom_types = eex['atom types'] = []
        masses = eex['masses'] = []

        for itype, xyz in zip(types.values(), coordinates):
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

    def eex_bond(self, eex):
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

    def eex_angle(self, eex):
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

    def eex_torsion(self, eex):
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

    def eex_out_of_plane(self, eex):
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

    def eex_bond_bond(self, eex):
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

    def eex_bond_angle(self, eex):
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

    def eex_torsion_middle_bond(self, eex):
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

    def eex_torsion_end_bond(self, eex):
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

    def eex_torsion_angle(self, eex):
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

    def eex_angle_torsion_angle(self, eex):
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

    def eex__3_bond_bond(self, eex):
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

    def eex_angle_angle(self, eex):
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

    def bond_parameters(self, i, j):
        """Return the bond parameters given two atoms types i and j

        Handle equivalences and automatic equivalences.
        """
        forms = self.atomtyping_engine.forcefield.ff['terms']['bond']
        # parameter directly available
        for form in forms:
            key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_bond', (i, j))
            if key in self.atomtyping_engine.forcefield.ff[form]:
                return ('explicit', key, form, self.atomtyping_engine.forcefield.ff[form][key])

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff:
            ieq = self.atomtyping_engine.forcefield.ff['equivalence'][i]['bond']
            jeq = self.atomtyping_engine.forcefield.ff['equivalence'][j]['bond']
            key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_bond', (ieq, jeq))
            for form in forms:
                if key in self.atomtyping_engine.forcefield.ff[form]:
                    return ('equivalent', key, form, self.atomtyping_engine.forcefield.ff[form][key])

        # try automatic equivalences
        if 'auto_equivalence' in self.atomtyping_engine.forcefield.ff:
            iauto = self.atomtyping_engine.forcefield.ff['auto_equivalence'][i]['bond']
            jauto = self.atomtyping_engine.forcefield.ff['auto_equivalence'][j]['bond']
            key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_bond', (iauto, jauto))
            for form in forms:
                if key in self.atomtyping_engine.forcefield.ff[form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff[form][key])
        import pdb; pdb.set_trace()
        raise RuntimeError('No bond parameters for {}-{}'.format(i, j))

    def angle_parameters(self, i, j, k):
        """Return the angle parameters given three atom types

        Handle equivalences and automatic equivalences.
        """

        forms = self.atomtyping_engine.forcefield.ff['terms']['angle']

        for form in forms:
            # parameters directly available
            result = self._angle_parameters_helper(i, j, k, self.atomtyping_engine.forcefield.ff[form])
            if result is not None:
                return ('explicit', result[0], form, result[2])

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff:
            ieq = self.atomtyping_engine.forcefield.ff['equivalence'][i]['angle']
            jeq = self.atomtyping_engine.forcefield.ff['equivalence'][j]['angle']
            keq = self.atomtyping_engine.forcefield.ff['equivalence'][k]['angle']
            for form in forms:
                result = self._angle_parameters_helper(
                    ieq, jeq, keq, self.atomtyping_engine.forcefield.ff[form]
                )
                if result is not None:
                    return ('equivalent', result[0], form, result[2])

        # try automatic equivalences
        if 'auto_equivalence' in self.atomtyping_engine.forcefield.ff:
            iauto = self.atomtyping_engine.forcefield.ff['auto_equivalence'][i]['angle_end_atom']
            jauto = self.atomtyping_engine.forcefield.ff['auto_equivalence'][j]['angle_center_atom']
            kauto = self.atomtyping_engine.forcefield.ff['auto_equivalence'][k]['angle_end_atom']
            key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                'like_angle', (iauto, jauto, kauto)
            )
            for form in forms:
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

            # try wildcards, which may have numerical precidence
            # Find all the single-sided wildcards, realizing that the
            # triplet might be flipped.
            for form in forms:
                left = []
                right = []
                for key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    if key[0] == '*' or key[2] == '*':
                        continue
                    if jauto == key[1]:
                        if kauto == key[2] and key[0][0] == '*':
                            left.append(key[0])
                        if kauto == key[0] and key[2][0] == '*':
                            left.append(key[2])
                        if iauto == key[0] and key[2][0] == '*':
                            right.append(key[2])
                        if iauto == key[2] and key[0][0] == '*':
                            right.append(key[0])
                if len(left) > 0:
                    if len(right) == 0:
                        key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                            'like_angle', (left[0], jauto, kauto)
                        )
                        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                            return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])
                    else:
                        if left[0] < right[0]:
                            key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                                'like_angle', (left[0], jauto, kauto)
                            )
                            if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                                return (
                                    'automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key]
                                )
                        else:
                            key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                                'like_angle', (iauto, jauto, right[0])
                            )
                            if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                                return (
                                    'automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key]
                                )
                elif len(right) > 0:
                    key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                        'like_angle', (iauto, jauto, right[0])
                    )
                    if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                        return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

                key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                    'like_angle', ('*', jauto, kauto)
                )
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])
                key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                    'like_angle', (iauto, jauto, '*')
                )
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])
                key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                    'like_angle', ('*', jauto, '*')
                )
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

        raise RuntimeError('No angle parameters for {}-{}-{}'.format(i, j, k))

    def torsion_parameters(self, i, j, k, l):  # noqa: E741
        """Return the torsion parameters given four atoms types

        Handles equivalences and automatic equivalences and wildcards,
        with numerical precedences
        """

        forms = self.atomtyping_engine.forcefield.ff['terms']['torsion']

        # parameter directly available
        for form in forms:
            result = self._torsion_parameters_helper(i, j, k, l, self.atomtyping_engine.forcefield.ff[form])
            if result is not None:
                return ('explicit', result[0], form, result[2])

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['torsion']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['torsion']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['torsion']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['torsion']
            for form in forms:
                result = self._torsion_parameters_helper(
                    ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff[form]
                )
                if result is not None:
                    return ('equivalent', result[0], form, result[2])

        # try automatic equivalences
        if 'auto_equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            iauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][i]['torsion_end_atom']
            jauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][j]['torsion_center_atom']
            kauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][k]['torsion_center_atom']
            lauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][l]['torsion_end_atom']
            key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                'like_torsion', (iauto, jauto, kauto, lauto)
            )
            for form in forms:
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

                # try wildcards, which may have numerical precidence
                # Find all the single-sided wildcards, realizing that the
                # triplet might be flipped.
                left = []
                right = []
                for key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    if key[0] == '*' or key[3] == '*':
                        continue
                    if jauto == key[1] and kauto == key[2]:
                        if lauto == key[3] and key[0][0] == '*':
                            left.append(key[0])
                        if lauto == key[0] and key[3][0] == '*':
                            left.append(key[3])
                        if iauto == key[0] and key[3][0] == '*':
                            right.append(key[3])
                        if iauto == key[3] and key[0][0] == '*':
                            right.append(key[0])
                if len(left) > 0:
                    if len(right) == 0:
                        key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                            'like_torsion', (left[0], jauto, kauto, lauto)
                        )
                        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                            return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])
                    else:
                        if left[0] < right[0]:
                            key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                                'like_torsion', (left[0], jauto, kauto, lauto)
                            )
                            if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                                return (
                                    'automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key]
                                )
                        else:
                            key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                                'like_torsion',
                                (iauto, jauto, kauto, right[0])
                            )
                            if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                                return (
                                    'automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key]
                                )
                elif len(right) > 0:
                    key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                        'like_torsion', (iauto, jauto, kauto, right[0])
                    )
                    if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                        return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

                key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                    'like_torsion', (iauto, jauto, kauto, '*')
                )
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])
                key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                    'like_torsion', ('*', jauto, kauto, lauto)
                )
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])
                key, flipped = self.atomtyping_engine.forcefield.make_canonical(
                    'like_torsion', ('*', jauto, kauto, '*')
                )
                if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                    return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

        raise RuntimeError(
            'No torsion parameters for {}-{}-{}-{}'.format(i, j, k, l)
        )

    def _torsion_parameters_helper(self, i, j, k, l, section):  # noqa: E741
        """Return the torsion parameters given four atom types
        """

        # parameter directly available
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_torsion', (i, j, k, l))
        if key in section:
            return (key, flipped, section[key])

        # try wildcards
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_torsion', ('*', j, k, l))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_torsion', (i, j, k, '*'))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_torsion', ('*', j, k, '*'))
        if key in section:
            return (key, flipped, section[key])

        return None

    def oop_parameters(self, i, j, k, l, zero=False):  # noqa: E741
        """Return the oop parameters given four atoms types

        Handles equivalences and automatic equivalences and wildcards,
        with numerical precedences
        """

        forms = self.atomtyping_engine.forcefield.ff['terms']['out-of-plane']

        for form in forms:
            result = self._oop_parameters_helper(i, j, k, l, form)
            if result is not None:
                return ('explicit', result[0], form, result[1])

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['oop']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['oop']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['oop']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['oop']
            for form in forms:
                result = self._oop_parameters_helper(ieq, jeq, keq, leq, form)
                if result is not None:
                    return ('equivalent', result[0], form, result[1])

        # try automatic equivalences
        if 'auto_equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            iauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][i]['oop_end_atom']
            jauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][j]['oop_center_atom']
            kauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][k]['oop_end_atom']
            lauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][l]['oop_end_atom']
            for form in forms:
                result = self._oop_parameters_helper(
                    iauto, jauto, kauto, lauto, form
                )
                if result is not None:
                    return ('automatic', result[0], form, result[1])

        if zero:
            parameters = {'K': 0.0, 'Chi0': 0.0}
            return (
                'zeroed', ('*', '*', '*', '*'), 'wilson_out_of_plane',
                parameters
            )
        else:
            raise RuntimeError(
                'No out-of-plane parameters for {}-{}-{}-{}'.format(
                    i, j, k, l
                )
            )

    def _oop_parameters_helper(self, i, j, k, l, form):  # noqa: E741
        """Return the oop parameters given four atoms types

        Handles equivalences and automatic equivalences and wildcards,
        with numerical precedences
        """

        # parameter directly available
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', (i, j, k, l))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])

        # try wildcards
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', ('*', j, k, l))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', (i, j, '*', l))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', (i, j, k, '*'))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', ('*', j, '*', l))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', ('*', j, k, '*'))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', (i, j, '*', '*'))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_oop', ('*', j, '*', '*'))
        if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
            return (key, self.atomtyping_engine.forcefield.ff['terms'][form][key])

        return None

    def nonbond_parameters(self, i, j=None, form='nonbond(12-6)'):
        """Return the nondbond parameters given one or two atoms types i and j

        Handle equivalences
        """

        # parameter directly available
        if j is None:
            key = (i,)
        else:
            key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_bond', (i, j))
        if key in self.atomtyping_engine.forcefield.ff[form]:
            return ('explicit', key, form, self.atomtyping_engine.forcefield.ff[form][key])

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['nonbond']
            if j is None:
                key = (ieq,)
            else:
                jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['nonbond']
                key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_bond', (ieq, jeq))
            if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                return ('equivalent', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

        # try automatic equivalences
        if 'auto_equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            iauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][i]['nonbond']
            if j is None:
                key = (iauto,)
            else:
                jauto = self.atomtyping_engine.forcefield.ff['terms']['auto_equivalence'][j]['nonbond']
                key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_bond', (iauto, jauto))
            if key in self.atomtyping_engine.forcefield.ff['terms'][form]:
                return ('automatic', key, form, self.atomtyping_engine.forcefield.ff['terms'][form][key])

        if j is None:
            raise RuntimeError('No nonbond parameters for {}'.format(i))
        else:
            raise RuntimeError('No nonbond parameters for {}-{}'.format(i, j))

    def bond_bond_parameters(self, i, j, k, zero=False):
        """Return the bond-bond parameters given three atoms types

        Handle equivalences, and if zero=True, return zero valued
        parameters rather than raise an error
        """

        # Get the reference bond lengths...
        b1_type, b1_types, b1_form, b1_parameters = \
            self.bond_parameters(i, j)
        b2_type, b2_types, b2_form, b2_parameters = \
            self.bond_parameters(j, k)
        values = {'R10': b1_parameters['R0'], 'R20': b2_parameters['R0']}

        # parameters directly available
        result = self._angle_parameters_helper(i, j, k, self.atomtyping_engine.forcefield.ff['bond-bond'])
        if result is not None:
            if result[1]:
                values = {
                    'R10': b2_parameters['R0'],
                    'R20': b1_parameters['R0']
                }
            values.update(result[2])
            return ('explicit', result[0], 'bond-bond', values)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['angle']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['angle']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['angle']
            result = self._angle_parameters_helper(
                ieq, jeq, keq, self.atomtyping_engine.forcefield.ff['bond-bond']
            )
            if result is not None:
                if result[1]:
                    values = {
                        'R10': b2_parameters['R0'],
                        'R20': b1_parameters['R0']
                    }
                values.update(result[2])
                return ('equivalent', result[0], 'bond-bond', values)

        if zero:
            return (
                'zeroed', ('*', '*', '*'), 'bond-bond', {
                    'K': '0.0',
                    'R10': '1.5',
                    'R20': '1.5'
                }
            )
        else:
            raise RuntimeError(
                'No bond-bond parameters for {}-{}-{}'.format(i, j, k)
            )

    def _angle_parameters_helper(self, i, j, k, section):
        """Return the angle-like parameters given three atom types
        """

        # parameter directly available
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle', (i, j, k))
        if key in section:
            return (key, flipped, section[key])

        # try wildcards
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle', ('*', j, k))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle', (i, j, '*'))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle', ('*', j, '*'))
        if key in section:
            return (key, flipped, section[key])

        return None

    def bond_bond_1_3_parameters(self, i, j, k, l, zero=False):  # noqa: E741
        """Return the bond-bond_1_3 parameters given four atoms types

        Handles equivalences wildcards
        """
        # Get the reference bond lengths...
        b1_type, b1_types, b1_form, b1_parameters = \
            self.bond_parameters(i, j)
        b3_type, b3_types, b3_form, b3_parameters = \
            self.bond_parameters(k, l)
        values = {'R10': b1_parameters['R0'], 'R30': b3_parameters['R0']}

        # parameter directly available
        result = self._torsion_parameters_helper(
            i, j, k, l, self.atomtyping_engine.forcefield.ff['bond-bond_1_3']
        )
        if result is not None:
            if result[1]:
                values = {
                    'R10': b3_parameters['R0'],
                    'R30': b1_parameters['R0']
                }
            values.update(result[2])
            return ('explicit', result[0], 'bond-bond_1_3', values)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['torsion']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['torsion']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['torsion']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff['bond-bond_1_3']
            )
            if result is not None:
                if result[1]:
                    values = {
                        'R10': b3_parameters['R0'],
                        'R30': b1_parameters['R0']
                    }
                values.update(result[2])
                return ('equivalent', result[0], 'bond-bond_1_3', values)

        if zero:
            parameters = {'K': '0.0', 'R10': '1.5', 'R30': '1.5'}
            return (
                'equivalent', ('*', '*', '*', '*'), 'bond-bond_1_3', parameters
            )
        else:
            raise RuntimeError(
                'No bond-bond_1_3 parameters for ' +
                '{}-{}-{}-{}'.format(i, j, k, l)
            )

    def bond_angle_parameters(self, i, j, k, zero=False):
        """Return the bond-angle parameters given three atoms types

        Handle equivalences, and if zero=True, return zero valued
        parameters rather than raise an error
        """

        # Get the reference bond lengths...
        b1_type, b1_types, b1_form, b1_parameters = \
            self.bond_parameters(i, j)
        b2_type, b2_types, b2_form, b2_parameters = \
            self.bond_parameters(j, k)

        # parameters directly available
        result = self._angle_parameters_helper(i, j, k, self.atomtyping_engine.forcefield.ff['bond-angle'])
        if result is not None:
            if result[1]:
                parameters = {
                    'reference': result[2]['reference'],
                    'K12': result[2]['K23'],
                    'K23': result[2]['K12'],
                    'R10': b2_parameters['R0'],
                    'R20': b1_parameters['R0']
                }
                ii, jj, kk = result[0]
                return ('explicit', (kk, jj, ii), 'bond-angle', parameters)
            else:
                parameters = dict(**result[2])
                parameters['R10'] = b1_parameters['R0']
                parameters['R20'] = b2_parameters['R0']
                return ('explicit', result[0], 'bond-angle', parameters)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['angle']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['angle']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['angle']
            result = self._angle_parameters_helper(
                ieq, jeq, keq, self.atomtyping_engine.forcefield.ff['bond-angle']
            )
            if result is not None:
                if result[1]:
                    parameters = {
                        'reference': result[2]['reference'],
                        'K12': result[2]['K23'],
                        'K23': result[2]['K12'],
                        'R10': b2_parameters['R0'],
                        'R20': b1_parameters['R0']
                    }
                    ii, jj, kk = result[0]
                    return (
                        'equivalent', (kk, jj, ii), 'bond-angle', parameters
                    )
                else:
                    parameters = dict(**result[2])
                    parameters['R10'] = b1_parameters['R0']
                    parameters['R20'] = b2_parameters['R0']
                    return ('equivalent', result[0], 'bond-angle', parameters)

        if zero:
            return (
                'zeroed', ('*', '*', '*'), 'bond-angle', {
                    'K12': '0.0',
                    'K23': '0.0',
                    'R10': '1.5',
                    'R20': '1.5'
                }
            )
        else:
            raise RuntimeError(
                'No bond-angle parameters for {}-{}-{}'.format(i, j, k)
            )

    def angle_angle_parameters(self, i, j, k, l, zero=False):  # noqa: E741
        """Return the angle_angle parameters given four atoms types

        Handles equivalences and wildcards
        """
        # Get the reference bond angles...
        a1_type, a1_types, a1_form, a1_parameters = \
            self.angle_parameters(i, j, k)
        a2_type, a2_types, a2_form, a2_parameters = \
            self.angle_parameters(k, j, l)
        Theta10 = a1_parameters['Theta0']
        Theta20 = a2_parameters['Theta0']
        values = {'Theta10': Theta10, 'Theta20': Theta20}

        # parameter directly available
        result = self._angle_angle_parameters_helper(
            i, j, k, l, self.atomtyping_engine.forcefield.ff['terms']['angle-angle']
        )
        if result is not None:
            if result[1]:
                values = {'Theta10': Theta20, 'Theta20': Theta10}
                values.update(result[2])
                ii, jj, kk, ll = result[0]
                return ('explicit', (ll, jj, kk, ii), 'angle-angle', values)
            else:
                values.update(result[2])
                return ('explicit', result[0], 'angle-angle', values)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['angle']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['angle']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['angle']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['angle']
            result = self._angle_angle_parameters_helper(
                ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff['terms']['angle-angle']
            )
            if result is not None:
                if result[1]:
                    values = {'Theta10': Theta20, 'Theta20': Theta10}
                    values.update(result[2])
                    ii, jj, kk, ll = result[0]
                    return (
                        'equivalent', (ll, jj, kk, ii), 'angle-angle', values
                    )
                else:
                    values.update(result[2])
                    return ('equivalent', result[0], 'angle-angle', values)

        if zero:
            parameters = {'K': 0.0, 'Theta10': '109.0', 'Theta20': '109.0'}
            return ('zeroed', ('*', '*', '*', '*'), 'angle-angle', parameters)
        else:
            raise RuntimeError(
                'No angle-angle parameters for {}-{}-{}-{}'.format(i, j, k, l)
            )

    def _angle_angle_parameters_helper(
        self,
        i,
        j,
        k,
        l,  # noqa: E741
        section
    ):
        """Return the torsion parameters given four atom types
        """

        # parameter directly available
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle-angle', (i, j, k, l))
        if key in section:
            return (key, flipped, section[key])

        # try wildcards
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle-angle', ('*', j, k, l))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical('like_angle-angle', (i, j, k, '*'))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical(
            'like_angle-angle', ('*', j, k, '*')
        )
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.atomtyping_engine.forcefield.make_canonical(
            'like_angle-angle', ('*', j, '*', '*')
        )
        if key in section:
            return (key, flipped, section[key])

        return None

    def end_bond_torsion_3_parameters(
        self,
        i,
        j,
        k,
        l,  # noqa: E741
        zero=False
    ):
        """Return the end bond - torsion_3 parameters given four atom types

        Handle equivalences
        """
        # Get the reference bond lengths...
        b1_type, b1_types, b1_form, b1_parameters = \
            self.bond_parameters(i, j)
        b2_type, b2_types, b2_form, b2_parameters = \
            self.bond_parameters(k, l)
        values = {'R0_L': b1_parameters['R0'], 'R0_R': b2_parameters['R0']}

        # parameters directly available
        result = self._torsion_parameters_helper(
            i, j, k, l, self.atomtyping_engine.forcefield.ff['terms']['end_bond-torsion_3']
        )
        if result is not None:
            if result[1]:
                parameters = {
                    'reference': result[2]['reference'],
                    'V1_L': result[2]['V1_R'],
                    'V2_L': result[2]['V2_R'],
                    'V3_L': result[2]['V3_R'],
                    'V1_R': result[2]['V1_L'],
                    'V2_R': result[2]['V2_L'],
                    'V3_R': result[2]['V3_L'],
                    'R0_L': b2_parameters['R0'],
                    'R0_R': b1_parameters['R0']
                }
                ii, jj, kk, ll = result[0]
                return (
                    'explicit', (ll, kk, jj, ii), 'end_bond-torsion_3',
                    parameters
                )
            else:
                parameters = dict(**result[2])
                parameters.update(values)
                return (
                    'explicit', result[0], 'end_bond-torsion_3', parameters
                )

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['torsion']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['torsion']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['torsion']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff['terms']['end_bond-torsion_3']
            )
            if result is not None:
                if result[1]:
                    parameters = {
                        'reference': result[2]['reference'],
                        'V1_L': result[2]['V1_R'],
                        'V2_L': result[2]['V2_R'],
                        'V3_L': result[2]['V3_R'],
                        'V1_R': result[2]['V1_L'],
                        'V2_R': result[2]['V2_L'],
                        'V3_R': result[2]['V3_L'],
                        'R0_L': b2_parameters['R0'],
                        'R0_R': b1_parameters['R0']
                    }
                    ii, jj, kk, ll = result[0]
                    return (
                        'equivalent', (ll, kk, jj, ii), 'end_bond-torsion_3',
                        parameters
                    )
                else:
                    parameters = dict(**result[2])
                    parameters.update(values)
                    return (
                        'equivalent', result[0], 'end_bond-torsion_3',
                        parameters
                    )

        if zero:
            parameters = {
                'V1_L': '0.0',
                'V2_L': '0.0',
                'V3_L': '0.0',
                'V1_R': '0.0',
                'V2_R': '0.0',
                'V3_R': '0.0',
                'R0_L': '1.5',
                'R0_R': '1.5'
            }
            return (
                'zeroed', ('*', '*', '*', '*'), 'end_bond-torsion_3',
                parameters
            )
        else:
            raise RuntimeError(
                'No end_bond-torsion_3 parameters for ' +
                '{}-{}-{}-{}'.format(i, j, k, l)
            )

    def middle_bond_torsion_3_parameters(
        self,
        i,
        j,
        k,
        l,  # noqa: E741
        zero=False
    ):
        """Return the middle bond - torsion_3 parameters given four atom types

        Handle equivalences
        """
        # Get the reference bond lengths...
        b1_type, b1_types, b1_form, b1_parameters = \
            self.bond_parameters(j, k)
        values = {'R0': b1_parameters['R0']}

        # parameters directly available
        result = self._torsion_parameters_helper(
            i, j, k, l, self.atomtyping_engine.forcefield.ff['terms']['middle_bond-torsion_3']
        )
        if result is not None:
            values.update(result[2])
            return ('explicit', result[0], 'middle_bond-torsion_3', values)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['torsion']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['torsion']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['torsion']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff['terms']['middle_bond-torsion_3']
            )
            if result is not None:
                values.update(result[2])
                return (
                    'equivalent', result[0], 'middle_bond-torsion_3', values
                )

        if zero:
            return (
                'zeroed', ('*', '*', '*', '*'), 'middle_bond-torsion_3', {
                    'R0': '1.5',
                    'V1': '0.0',
                    'V2': '0.0',
                    'V3': '0.0'
                }
            )
        else:
            raise RuntimeError(
                'No middle_bond-torsion_3 parameters for ' +
                '{}-{}-{}-{}'.format(i, j, k, l)
            )

    def angle_torsion_3_parameters(self, i, j, k, l, zero=False):  # noqa: E741
        """Return the angle - torsion_3 parameters given four atom types

        Handle equivalences
        """
        # Get the reference bond angles...
        a1_type, a1_types, a1_form, a1_parameters = \
            self.angle_parameters(i, j, k)
        a2_type, a2_types, a2_form, a2_parameters = \
            self.angle_parameters(j, k, l)
        values = {
            'Theta0_L': a1_parameters['Theta0'],
            'Theta0_R': a2_parameters['Theta0']
        }

        # parameters directly available
        result = self._torsion_parameters_helper(
            i, j, k, l, self.atomtyping_engine.forcefield.ff['terms']['angle-torsion_3']
        )
        if result is not None:
            if result[1]:
                parameters = {
                    'reference': result[2]['reference'],
                    'V1_L': result[2]['V1_R'],
                    'V2_L': result[2]['V2_R'],
                    'V3_L': result[2]['V3_R'],
                    'V1_R': result[2]['V1_L'],
                    'V2_R': result[2]['V2_L'],
                    'V3_R': result[2]['V3_L'],
                    'Theta0_L': a2_parameters['Theta0'],
                    'Theta0_R': a1_parameters['Theta0']
                }
                ii, jj, kk, ll = result[0]
                return (
                    'explicit', (ll, kk, jj, ii), 'angle-torsion_3', parameters
                )
            else:
                parameters = dict(**result[2])
                parameters.update(values)
                return ('explicit', result[0], 'angle-torsion_3', parameters)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['torsion']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['torsion']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['torsion']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff['terms']['angle-torsion_3']
            )
            if result is not None:
                if result[1]:
                    parameters = {
                        'reference': result[2]['reference'],
                        'V1_L': result[2]['V1_R'],
                        'V2_L': result[2]['V2_R'],
                        'V3_L': result[2]['V3_R'],
                        'V1_R': result[2]['V1_L'],
                        'V2_R': result[2]['V2_L'],
                        'V3_R': result[2]['V3_L'],
                        'Theta0_L': a2_parameters['Theta0'],
                        'Theta0_R': a1_parameters['Theta0']
                    }
                    ii, jj, kk, ll = result[0]
                    return (
                        'equivalent', (ll, kk, jj, ii), 'angle-torsion_3',
                        parameters
                    )
                else:
                    parameters = dict(**result[2])
                    parameters.update(values)
                    return (
                        'equivalent', result[0], 'angle-torsion_3', parameters
                    )

        if zero:
            parameters = {
                'V1_L': '0.0',
                'V2_L': '0.0',
                'V3_L': '0.0',
                'V1_R': '0.0',
                'V2_R': '0.0',
                'V3_R': '0.0',
                'Theta0_L': '109.0',
                'Theta0_R': '109.0'
            }
            return (
                'zeroed', ('*', '*', '*', '*'), 'angle-torsion_3', parameters
            )
        else:
            raise RuntimeError(
                'No angle-torsion_3 parameters for ' +
                '{}-{}-{}-{}'.format(i, j, k, l)
            )

    def angle_angle_torsion_1_parameters(
        self,
        i,
        j,
        k,
        l,  # noqa: E741
        zero=False
    ):
        """Return the angle - angle - torsion_1 parameters given four atom types

        Handle equivalences
        """
        # Get the reference bond angles...
        a1_type, a1_types, a1_form, a1_parameters = \
            self.angle_parameters(i, j, k)
        a2_type, a2_types, a2_form, a2_parameters = \
            self.angle_parameters(j, k, l)
        values = {
            'Theta0_L': a1_parameters['Theta0'],
            'Theta0_R': a2_parameters['Theta0']
        }

        # parameters directly available
        result = self._torsion_parameters_helper(
            i, j, k, l, self.atomtyping_engine.forcefield.ff['terms']['angle-angle-torsion_1']
        )
        if result is not None:
            values.update(result[2])
            return ('explicit', result[0], 'angle-angle-torsion_1', values)

        # try equivalences
        if 'equivalence' in self.atomtyping_engine.forcefield.ff['terms']:
            ieq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][i]['torsion']
            jeq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][j]['torsion']
            keq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][k]['torsion']
            leq = self.atomtyping_engine.forcefield.ff['terms']['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.atomtyping_engine.forcefield.ff['terms']['angle-angle-torsion_1']
            )
            if result is not None:
                values.update(result[2])
                return (
                    'equivalent', result[0], 'angle-angle-torsion_1', values
                )

        if zero:
            parameters = {'Theta0_L': '109.0', 'Theta0_R': '109.0', 'K': '0.0'}
            return (
                'zeroed', ('*', '*', '*', '*'), 'angle-angle-torsion_1',
                parameters
            )
        else:
            raise RuntimeError(
                'No angle-angle-torsion_1 parameters for ' +
                '{}-{}-{}-{}'.format(i, j, k, l)
            )

    def eex_charges(self, eex):
        """Do nothing routine since charges are handled by the increments."""
        pass

    def eex_increment(self, eex):
        """Get the charges for the structure

        If they do not exists on the structure, they are created
        using the bond increments and saved on the structure"""

        logger.debug('entering eex_increment')

        key = f'charges_{self.atomtyping_engine.forcefield.name}'
        if key in self.configuration.atoms:
            eex['charges'] = [*self.configuration.atoms[key]]
        else:
            raise RuntimeError('No charges on system!')

        logger.debug('leaving eex_increment')


    def eex_pair(self, eex):
        """Create the pair (non-bond) portion of the energy expression"""
        logger.debug('In eex_pair')
        types = self.topology['types']

        found = False
        nonbond_types = ('nonbond(12-6)', 'nonbond(9-6)') 

        for k, v in self.atomtyping_engine.forcefield.ff['terms'].items():
            intersection = set(nonbond_types) & set(v) 
            if len(intersection) > 0:
                pair_type = intersection.pop()
                found = True
                break

        if not found:
            raise RuntimeError('Error finding pair_type in eex_pair')

        result = eex['nonbonds'] = []
        parameters = eex['nonbond parameters'] = []
        for k, itype in types.items():
            parameters_type, real_types, form, parameter_values = \
                self.nonbond_parameters(itype, form=pair_type)
            new_value = (
                form, parameter_values, (itype,), parameters_type, real_types
            )
            index = None
            for value, count in zip(parameters, range(1, len(parameters) + 1)):
                if new_value == value:
                    index = count
                    break
            if index is None:
                parameters.append(new_value)
                index = len(parameters)
            result.append(index)
        eex['n_nonbonds'] = len(result)
        eex['n_nonbond_types'] = len(parameters)
