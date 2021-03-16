# -*- coding: utf-8 -*-

# Don't we need some versioning function in the sort? (about line 1200)

"""Main class for handling forcefields"""

from enum import Enum
import json
import logging
import os.path
import packaging.version
import pprint
from .forcefield_metadata import metadata
import seamm_util
from seamm_util import Q_

logger = logging.getLogger(__name__)


class Forcefield(object):

    def __init__(self, filename=None, version=None):
        """
        Read, write, and use a forcefield from a biosym file

        The Forcefield object is the main interface for working with
        forcefields. It provides methods to read and write
        forcefields, to assign the forcefield to a molecule, as well
        as to get parameters for bonds, angles, etc.

        Args:
            filename ('str', optional): An optional filename for the forcefield
            fftype ('str', optional): An optional type for the
                forcefield. If not given and a forcefield is read, the
                code will try to divine the type of forcefield.
        """

        self.keep_lines = False

        self.name = None

        self._filename = None

        self.data = {}

        self.ff = {}

        self.filename = filename

        self.version = None

        if version is not None:
            self.version = packaging.version.Version(version)

        fforms = self.data['forcefield'][self.name]['parameters']
        self.ff['terms'] = []
        self.ff['modifiers'] = {}
        self.ff['functional_forms'] = {}

        for fform in fforms:

            versions = sorted(fforms[fform].keys(), reverse=True)

            if version is None:
                key = versions[0]
            else:
                key = None
                for value in versions:
                    if value <= self.version:
                        key = value
                        break
                if key is None:
                    raise RuntimeError(
                        "Cannot find version '{}'".format(self.version) +
                        " for functional form '{}'".format(fform) +
                        " of forcefield '{}'".format(forcefield)
                    )

            self.ff['functional_forms'][fform] = fforms[fform][key]

            if fform in metadata:

                term = metadata[fform]['topology']['type']

                if term in self.ff['terms']:
                    self.ff['terms'].append(fform)
                else:
                    self.ff['terms'].append(fform)

            for fform in self.ff['functional_forms']:
                self._get_parameters(fform, self.version)

    @property
    def filename(self):
        """'str' name of file for this forcefield.

        When the filename is set, if the file exists it is read. If it
        does not exist, it is created and initialized as a forcefield
        file. The type of the file may be given by self.fftype; if
        not, the code tries the divine the type of the forcefield. The
        default type for new forcefields is the Biosym .frc format.

        If the filename is changed the object is reset.
        """
        return self._filename

    @filename.setter
    def filename(self, value):
        if not value:
            self.clear()
            self._filename = None
        else:
            if value == self._filename:
                return

            if os.path.isfile(value):
                self.clear()

                self.name = None 
                self._filename = value

                with seamm_util.Open(self._filename, 'r') as fd:
                    self._read_biosym_ff(fd)
            else:
                self._filename = value
                self._create_file()

    def clear(self):
        """
        Reset the object to its initial, empty, state
        """
        self._filename = None
        self.data = {}
        self.ff = {}
        self.name = None

    def _read_biosym_ff(self, fd):
        """
        Read and parse a forcefield in Biosym's format

        Args:
            fd (file object): the file handle
        """


        self.name = None
        self.data = {
            'forcefield': {},
        }

        try:
            # Read and process the first line, which should say
            # what the file is e.g. '!BIOSYM forcefield 1'
            line = next(fd)
            if line[0] == '!' and len(line.split()) == 3:
                file_variant, file_type, version = line[1:].split()
                logger.info(
                    "reading '{}', a {} file from {}, version {}".format(
                        self.filename, file_type, file_variant, version
                    )
                )
            else:
                logger.warning(
                    "reading '{}', expected a header line but got\n\t'{}'"
                    .format(self.filename, line)
                )

            # Read the rest of the file, processing the '#'
            # delimited sections
            for line in fd:
                line = line.strip()

                # Empty and comment lines
                if line == '' or line[0] == '!':
                    continue

                if line[0] == '#':
                    # fd.push()
                    words = line[1:].split()
                    section = words[0]

                    # Just ignore #end sections, as they simply close a section
                    if section == 'end' or section == 'version':
                        continue

                    if len(words) < 2:
                        logger.warning(
                            section + ' section does not have a label!\n\t' +
                            '\n\t'.join(fd.stack())
                        )
                        label = 'missing'
                        priority = 0
                    else:
                        label = words[1]
                        if len(words) > 2:
                            priority = float(words[2])
                        else:
                            priority = 0

                    logger.debug('reading ff section ' + section)
                    
                    result = self._read_biosym_section(fd)

                    result['section'] = section
                    result['label'] = label
                    result['priority'] = priority

                    # Parse the data, looking for specialized implementations
                    if 'nonbond' in section:
                        method = '_parse_biosym_nonbonds'
                    else:
                        method = '_parse_biosym_' + section
                    logger.info(
                        "Parsing forcefield section '" + section + "'."
                    )
                    if method in Forcefield.__dict__:
                        Forcefield.__dict__[method](self, result)
                    elif section in metadata:
                        self._parse_biosym_section(result)
                    else:
                        logger.warning('Cannot find parser for ' + section)

        except IOError:
            logger.exception(
                "Encountered I/O error opening '{}'".format(self.filename)
            )
            raise

    def _read_biosym_section(self, fd):
        """
        Read the body of a section of the forcefield

        Keeps tracks of comments ('!'), annotations ('>'), and modifiers ('@'),
        returning a dictionary with them plus tte raw lines of data
        """
        result = {
            'comments': [],
            'lines': [],
            'annotations': [],
            'modifiers': []
        }


        for line in fd:
            line = line.strip()

            # Empty and comment lines
            if line == '':
                continue

            if line[0] == '!':
                result['comments'].append(line[1:])
                continue

            if line[0] == '#':
                # At the end of the section, push the line back so the
                # main reader handles it and return the dict with the
                # data
                fd.push()
                return result

            if line[0] == '>':
                # An annotation
                result['annotations'].append(line[1:])
                continue

            if line[0] == '@':
                # A modifier such as units or form
                result['modifiers'].append(line[1:])
                continue

            # Must be a line of data! :-)
            result['lines'].append(line)
        return result 

    def _parse_biosym_define(self, data):
        """
        Process a forcefield definition section

        #define cff91

        !Ver Ref		Function	     Label
        !--- ---    ------------------------------   ------
         1.0  1     atom_types                       cff91
         1.0  1     equivalence                      cff91
        ...
        """
        section = 'forcefield'
        ff_name = data['label']
        self.name = ff_name

        if section not in self.data:
            self.data[section] = {}
        self.data[section][ff_name] = data
        sections = self.data[section][ff_name]['parameters'] = {}

        for line in data['lines']:
            words = line.split()
            if len(words) < 4:
                logger.error(
                    "In a define section for {}, the line is too short:"
                    .format(ff_name)
                )
                logger.error("    " + line)
            else:
                version, reference, functional_form = words[0:3]
                labels = words[3:]
                if functional_form not in sections:
                    sections[functional_form] = {}
                V = packaging.version.Version(version)
                sections[functional_form][V] = {
                    'version': version,
                    'reference': reference,
                    'sections': labels
                }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_atom_types(self, data):
        """
        Process the atom types

        #atom_types           cff91

        > Atom type definitions for any variant of cff91
        > Masses from CRC 1973/74 pages B-250.

        !Ver Ref  Type     Mass      Element   connection   Comment
        !--- ---  -----  ----------  -------   ----------   -------------------
        2.1 11   Ag     107.86800     Ag          0        Silver metal
        2.1 11   Al      26.98200     Al          0        Aluminium metal
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        atom_types = self.data[section][label]['parameters'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference, atom_type, mass, element, connections = words[
                0:6]
            comment = ' '.join(words[6:])
            if atom_type not in atom_types:
                atom_types[atom_type] = {}
            V = packaging.version.Version(version)
            if V in atom_types[atom_type]:
                msg = "atom type '{}' defined more than ".format(atom_type) + \
                      "once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            atom_types[atom_type][V] = {
                'reference': reference,
                'mass': mass,
                'element': element,
                'connections': connections,
                'comment': comment
            }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_equivalence(self, data):
        """
        Process the atom type equivalences

        #equivalence          cff91

        !                      Equivalences
        !       ------------------------------------------
        !Ver Ref  Type   NonB   Bond   Angle  Torsion  OOP
        !--- ---  -----  -----  -----  -----  -------  -----
        2.1 11   Ag     Ag     Ag     Ag     Ag       Ag
        2.1 11   Al     Al     Al     Al     Al       Al
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        equivalences = self.data[section][label]['parameters'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference, atom_type, nonbond, bond, angle, \
                torsion, oop = words
            if atom_type not in equivalences:
                equivalences[atom_type] = {}
            V = packaging.version.Version(version)
            if V in equivalences[atom_type]:
                msg = "atom type '{}' defined more than ".format(atom_type) + \
                      "once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            equivalences[atom_type][V] = {
                'reference': reference,
                'nonbond': nonbond,
                'bond': bond,
                'angle': angle,
                'torsion': torsion,
                'oop': oop
            }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_auto_equivalence(self, data):
        """
        Process the atom type equivalences for automatic types

        #auto_equivalence     cff91_auto

        !                      Equivalences
        !       ------------------------------------------
        !Ver  Ref   Type  NonB Bond   Bond     Angle    Angle     Torsion   Torsion      OOP      OOP
        !                       Inct           End atom Apex atom End Atoms Center Atoms End Atom Center Atom
        !---- ---   ----  ---- ------ ----  ---------- --------- --------- -----------  -------- -----------
        2.0  1     Br    Br   Br     Br_   Br_        Br_       Br_       Br_          Br_      Br_
        2.0  1     Cl    Cl   Cl     Cl_   Cl_        Cl_       Cl_       Cl_          Cl_      Cl_
        ...
        """  # noqa: E501
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        equivalences = self.data[section][label]['parameters'] = {}

        for line in data['lines']:
            words = line.split()
            version, reference, atom_type, nonbond, bond_increment, bond, \
                angle_end_atom, angle_center_atom, torsion_end_atom, \
                torsion_center_atom, oop_end_atom, oop_center_atom = words
            if atom_type not in equivalences:
                equivalences[atom_type] = {}
            V = packaging.version.Version(version)
            if V in equivalences[atom_type]:
                msg = "atom type '{}' defined more than ".format(atom_type) + \
                      "once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            equivalences[atom_type][V] = {
                'reference': reference,
                'nonbond': nonbond,
                'bond_increment': bond_increment,
                'bond': bond,
                'angle_end_atom': angle_end_atom,
                'angle_center_atom': angle_center_atom,
                'torsion_end_atom': torsion_end_atom,
                'torsion_center_atom': torsion_center_atom,
                'oop_end_atom': oop_end_atom,
                'oop_center_atom': oop_center_atom
            }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_bond_increments(self, data):
        """
        Process the bond increments

        #bond_increments      cff91_auto

        !Ver Ref    I     J     DeltaIJ   DeltaJI
        !--- ---  ----- -----   -------   -------
        2.1 11   Ag    Ag       0.0000   0.0000
        2.1 11   Al    Al       0.0000   0.0000
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)

        self.data[section][label] = data

        parameters = data['parameters'] = {}

        # Copy in the metadata about this functional form
        data.update(metadata[section])

        for line in data['lines']:
            words = line.split()
            version, reference, i, j, deltaij, deltaji = words
            # order canonically, i<j
            if i > j:
                i, j = j, i
                deltaij, deltaji = deltaji, deltaij
            key = (i, j)
            if key not in parameters:
                parameters[key] = {}
            V = packaging.version.Version(version)
            if V in parameters[key]:
                msg = "bond increment '{}' '{}' defined more ".format(i, j) + \
                      "than once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)
            parameters[key][V] = {
                'reference': reference,
                'deltaij': deltaij,
                'deltaji': deltaji
            }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_templates(self, data):
        """
        Process the templates, which are simply json

        #templates pcff
        "c": {
            "2017.12.15": {
                "smarts": [
                    "[CX4:1]"
                ],
                "description": "generic SP3 carbon",
                "overrides": []
            }
        },
        "c3": {
            "2017.12.15": {
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data

        data['parameters'] = json.loads('\n'.join(data['lines']))

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_reference(self, data):
        """
        Process a 'reference' section, which looks like

        #reference 1
        @Author Biosym Technologies inc
        @Date 25-December-91
        cff91 forcefield created
        December 1991

        """
        section = data['section']
        label = data['label']

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)
        self.data[section][label] = data
        data['reference'] = data['lines']

        if not self.keep_lines:
            del data['lines']

    def make_canonical(self, symmetry, atom_types):
        """
        Using the symmetry, order the atom_types canonically
        """

        n = len(atom_types)
        flipped = False
        if n == 1:
            i = atom_types[0]
            return ((i,), flipped)
        elif n == 2:
            i, j = atom_types[0:2]
            if symmetry == 'like_bond':
                # order canonically, i<j
                if i > j:
                    i, j = j, i
                    flipped = True
                return ((i, j), flipped)
        elif n == 3:
            i, j, k = atom_types[0:3]
            if symmetry == 'like_angle':
                # order canonically, i<k
                if i > k:
                    i, k = k, i
                    flipped = True
                return ((i, j, k), flipped)
        elif n == 4:
            i, j, k, l = atom_types[0:4]  # noqa: E741
            if symmetry == 'like_torsion':
                # order canonically, j<k; i<l if j==k
                if j == k and i > l:
                    i, l = l, i  # noqa: E741
                    flipped = True
                elif j > k:
                    i, j, k, l = l, k, j, i  # noqa: E741
                    flipped = True
                return ((i, j, k, l), flipped)
            elif symmetry == 'like_oop':
                # j is central atom
                # order canonically, i<k<l; i=k<l or i<k=l
                i, k, l = sorted((i, k, l))  # noqa: E741
                flipped = [i, j, k, l] != atom_types
                return ((i, j, k, l), flipped)
            elif symmetry == 'like_angle-angle':
                # order canonically, i<l;
                if i > l:
                    i, l = l, i  # noqa: E741
                    flipped = True
                return ((i, j, k, l), flipped)

    def _parse_biosym_nonbonds(self, data):
        """
        Process the nonbond parameters, accounting for different expressions
        and units.

        For example:
            #nonbond(12-6) spc

            > E = (A/r)^12 - (B/r)^6
            >
            > where    r(ij) is the distance between atoms i and j

            @type A/r-B/r
            @units A (kJ/mol)**(1/12)*nm
            @units B (kJ/mol)**(1/6)*nm
            @combination geometric

            !   Ver    Ref    I            A              B
            !--------- ---  -------   -------------  -----------
            2020.05.13   3  o_spc       0.3428         0.37122
            2020.05.13   3  h_spc       0.0            0.0

            #nonbond(12-6) nacl

            >   E = eps * [(rmin/r)^12 - (rmin/r)^6]
            >
            > where    r is the distance between atoms i and j

            @type eps-rmin
            @units eps kcal/mol
            @units rmin angstrom
            @combination geometric

            !   Ver    Ref    I           Rmin        Epsilon
            !--------- ---  -------   -------------  -----------
            2020.05.13   4  na+         2.7275         0.0469
            2020.05.13   4  k+          3.5275         0.0870
            2020.05.13   5  cl-         4.5400         0.1500


        """  # nopep8
        logger.debug('Entering _parse_biosym_nonbonds')

        section = data['section']
        label = data['label']

        logger.debug('parsing section ' + section + ' with nonbond parser')
        logger.debug('  data keys: ' + str(data.keys()))

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)

        self.data[section][label] = data

        # Copy in the metadata about this functional form
        data.update(metadata[section])

        topology = data['topology']

        out1, out2 = data['constants']
        out1_units = out1[1]  # angstrom, nm
        out2_units = out2[1]  # kcal/mol, kJ/mol
        parameter_1 = out1[0]
        parameter_2 = out2[0]

        out_form = NonbondForms(topology['form'])

        # Default for parameters read in
        in_form = out_form
        in1_units = 'angstrom'
        in2_units = 'kcal/mol'

        # And see if there are modifiers
        for item in data['modifiers']:
            modifier = item.split()
            what = modifier[0]
            if what == 'type':
                try:
                    in_form = NonbondForms(modifier[1])
                except ValueError:
                    raise ValueError(
                        "Unrecognized nonbond form '" + modifier[1] + "'"
                    )
            elif what == 'units':
                which = modifier[1]
                if which in ['sigma', 'rmin', 'A']:
                    in1_units = modifier[2]
                elif which in ['eps', 'B']:
                    in2_units = modifier[2]
                else:
                    raise ValueError(
                        "Unrecognized nonbond parameter '" + which + "'"
                    )
            elif what == 'combination':
                pass
            else:
                raise ValueError(
                    "Unrecognized nonbond modifier '" + what + "' Should be "
                    "one of 'type', 'units' or 'combination'."
                )

        # Now that we know what we are getting as parameters (in1, in2),
        # and what we want (p1, p2), make it so!
        transform = Forcefield.nonbond_transformation(
            in_form, in1_units, in2_units, out_form, out1_units, out2_units
        )

        parameters = data['parameters'] = {}

        for line in data['lines']:
            version, reference, i, p1, p2 = line.split()

            key = (i,)
            if key not in parameters:
                parameters[key] = {}
            V = packaging.version.Version(version)
            if V in parameters[key]:
                msg = "value for '" + "' '".join(key) + " defined more " + \
                      "than once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)

            v1, v2 = transform(float(p1), float(p2))

            parameters[key][V] = {
                'reference': reference,
                parameter_1: v1,
                parameter_2: v2
            }

        if not self.keep_lines:
            del data['lines']

    def _parse_biosym_section(self, data):
        """
        Process the 1-term torsion parameters

        #torsion_1            cff91_auto

        > E = Kphi * [ 1 + cos(n*Phi - Phi0) ]

        !Ver Ref    I     J     K     L       KPhi     n     Phi0
        !--- ---  ----- ----- ----- -----   --------  ---  ---------
        2.0  2   *     c'_   c'_   *         0.4500    2   180.0000
        2.0  2   *     c'_   c=_   *         0.4500    2   180.0000
        ...
        """  # nopep8
        section = data['section']
        label = data['label']

        logger.debug('parsing section ' + section + ' with generic parser')
        logger.debug('  data keys: ' + str(data.keys()))

        if section not in self.data:
            self.data[section] = {}
        if label in self.data[section]:
            msg = "'{}' already defined in section '{}'".format(label, section)
            logger.error(msg)
            raise RuntimeError(msg)

        self.data[section][label] = data

        # Copy in the metadata about this functional form
        data.update(metadata[section])

        parameters = data['parameters'] = {}

        for line in data['lines']:
            
            words = line.split()
            version, reference = words[0:2]
            symmetry = data['topology']['symmetry']
            n_atoms = data['topology']['n_atoms']
            key, flipped = self.make_canonical(symmetry, words[2:2 + n_atoms])

            V = packaging.version.Version(version)

            if key not in parameters:
                term_idx = 0
                parameters[key]= {}
                parameters[key][V] = {}
                parameters[key][V]['reference'] = reference
                parameters[key][V]['terms'] = {}

            if V not in parameters[key]:
                parameters[key][V] = {}
                parameters[key][V]['reference'] = reference
                parameters[key][V]['terms'] = {}
            
            values = words[2 + n_atoms:]

            if 'fill' in data['topology']:
                n = data['topology']['fill']
                if n > 0:
                    if len(values) < 2 * n:
                        values.extend(values[0:n])
            if flipped and 'flip' in data['topology']:
                n = data['topology']['flip']
                if n > 0:
                    first = values[0:n]
                    values = values[n:2 * n]
                    values.extend(first)

            if term_idx in parameters[key][V]['terms']:
                msg = "value for '" + "' '".join(key) + " defined more " + \
                      "than once in section '{}'!".format(section)
                logger.error(msg)
                raise RuntimeError(msg)

            val_dict = {const[0]: val for const, val in zip(data['constants'], values)}

            parameters[key][V]['terms'][term_idx] = val_dict

            term_idx += 1

        #    for constant, value in zip(data['constants'], values):
        #        params[constant[0]] = value
        if not self.keep_lines:
            del data['lines']


    def _get_parameters(self, functional_form, version):
        """Select the correct version parameters from the sections for
        this functional form"""

        logger.debug('_get_parameters, form = ' + functional_form)
        sections = self.ff['functional_forms'][functional_form]['sections']

        logger.debug('  sections = ' + str(sections))

        newdata = self.ff[functional_form] = {}
        modifiers = self.ff['modifiers'][functional_form] = {}

        for section in sections:
            data = self.data[functional_form][section]['parameters']

            modifiers[section] = \
                self.data[functional_form][section]['modifiers']

            for item in data:

                # Don't we need some versioning function in the sort?
                versions = sorted(data[item].keys(), reverse=True)

                if version is None:
                    key = versions[0]
                else:
                    key = None
                    for value in versions:
                        if value <= version:
                            key = value
                            break
                if key is not None:
                    newdata[item] = data[item][key]


    def charges(self, i):
        """Return the charge given an atom type i

        Handle equivalences.
        """

        if 'charges' in self.ff:
            # parameter directly available
            key = (i,)
            if key in self.ff['charges']:
                parameters = {}
                parameters.update(self.ff['charges'][key])
                return ('explicit', key, 'charges', parameters)

            # try equivalences
            if 'equivalence' in self.ff:
                ieq = self.ff['equivalence'][i]['nonbond']
                key = (ieq,)
                if key in self.ff['charges']:
                    parameters = {}
                    parameters.update(self.ff['charges'][key])
                    return ('equivalent', key, 'charges', parameters)

        # return the default of zero
        parameters = {'Q': 0.0}
        return ('default', ('*',), 'charges', parameters)

    def bond_increments(self, i, j):
        """Return the bond increments given two atoms types i and j

        Handle automatic equivalences.
        """

        # parameter directly available
        key, flipped = self.make_canonical('like_bond', (i, j))
        if key in self.ff['bond_increments']:
            parameters = {}
            parameters.update(self.ff['bond_increments'][key])
            if flipped:
                parameters['deltaij'], parameters['deltaji'] = \
                    parameters['deltaji'], parameters['deltaij']
            return ('explicit', key, 'bond_increments', parameters)

        # try automatic equivalences
        if 'auto_equivalence' in self.ff:
            iauto = self.ff['auto_equivalence'][i]['bond_increment']
            jauto = self.ff['auto_equivalence'][j]['bond_increment']
            key, flipped = self.make_canonical('like_bond', (iauto, jauto))
            if key in self.ff['bond_increments']:
                parameters = {}
                parameters.update(self.ff['bond_increments'][key])
                if flipped:
                    parameters['deltaij'], parameters['deltaji'] = \
                        parameters['deltaji'], parameters['deltaij']
                return ('automatic', key, 'bond_increments', parameters)

        raise RuntimeError('No bond increments for {}-{}'.format(i, j))

    def bond_parameters(self, i, j):
        """Return the bond parameters given two atoms types i and j

        Handle equivalences and automatic equivalences.
        """

        forms = self.ff['terms']['bond']

        # parameter directly available
        for form in forms:
            key, flipped = self.make_canonical('like_bond', (i, j))
            if key in self.ff[form]:
                return ('explicit', key, form, self.ff[form][key])

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['bond']
            jeq = self.ff['equivalence'][j]['bond']
            key, flipped = self.make_canonical('like_bond', (ieq, jeq))
            for form in forms:
                if key in self.ff[form]:
                    return ('equivalent', key, form, self.ff[form][key])

        # try automatic equivalences
        if 'auto_equivalence' in self.ff:
            iauto = self.ff['auto_equivalence'][i]['bond']
            jauto = self.ff['auto_equivalence'][j]['bond']
            key, flipped = self.make_canonical('like_bond', (iauto, jauto))
            for form in forms:
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])

        raise RuntimeError('No bond parameters for {}-{}'.format(i, j))

    def angle_parameters(self, i, j, k):
        """Return the angle parameters given three atom types

        Handle equivalences and automatic equivalences.
        """

        forms = self.ff['terms']['angle']

        for form in forms:
            # parameters directly available
            result = self._angle_parameters_helper(i, j, k, self.ff[form])
            if result is not None:
                return ('explicit', result[0], form, result[2])

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['angle']
            jeq = self.ff['equivalence'][j]['angle']
            keq = self.ff['equivalence'][k]['angle']
            for form in forms:
                result = self._angle_parameters_helper(
                    ieq, jeq, keq, self.ff[form]
                )
                if result is not None:
                    return ('equivalent', result[0], form, result[2])

        # try automatic equivalences
        if 'auto_equivalence' in self.ff:
            iauto = self.ff['auto_equivalence'][i]['angle_end_atom']
            jauto = self.ff['auto_equivalence'][j]['angle_center_atom']
            kauto = self.ff['auto_equivalence'][k]['angle_end_atom']
            key, flipped = self.make_canonical(
                'like_angle', (iauto, jauto, kauto)
            )
            for form in forms:
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])

            # try wildcards, which may have numerical precidence
            # Find all the single-sided wildcards, realizing that the
            # triplet might be flipped.
            for form in forms:
                left = []
                right = []
                for key in self.ff[form]:
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
                        key, flipped = self.make_canonical(
                            'like_angle', (left[0], jauto, kauto)
                        )
                        if key in self.ff[form]:
                            return ('automatic', key, form, self.ff[form][key])
                    else:
                        if left[0] < right[0]:
                            key, flipped = self.make_canonical(
                                'like_angle', (left[0], jauto, kauto)
                            )
                            if key in self.ff[form]:
                                return (
                                    'automatic', key, form, self.ff[form][key]
                                )
                        else:
                            key, flipped = self.make_canonical(
                                'like_angle', (iauto, jauto, right[0])
                            )
                            if key in self.ff[form]:
                                return (
                                    'automatic', key, form, self.ff[form][key]
                                )
                elif len(right) > 0:
                    key, flipped = self.make_canonical(
                        'like_angle', (iauto, jauto, right[0])
                    )
                    if key in self.ff[form]:
                        return ('automatic', key, form, self.ff[form][key])

                key, flipped = self.make_canonical(
                    'like_angle', ('*', jauto, kauto)
                )
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])
                key, flipped = self.make_canonical(
                    'like_angle', (iauto, jauto, '*')
                )
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])
                key, flipped = self.make_canonical(
                    'like_angle', ('*', jauto, '*')
                )
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])

        raise RuntimeError('No angle parameters for {}-{}-{}'.format(i, j, k))

    def torsion_parameters(self, i, j, k, l):  # noqa: E741
        """Return the torsion parameters given four atoms types

        Handles equivalences and automatic equivalences and wildcards,
        with numerical precedences
        """

        forms = self.ff['terms']['torsion']

        # parameter directly available
        for form in forms:
            result = self._torsion_parameters_helper(i, j, k, l, self.ff[form])
            if result is not None:
                return ('explicit', result[0], form, result[2])

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['torsion']
            jeq = self.ff['equivalence'][j]['torsion']
            keq = self.ff['equivalence'][k]['torsion']
            leq = self.ff['equivalence'][l]['torsion']
            for form in forms:
                result = self._torsion_parameters_helper(
                    ieq, jeq, keq, leq, self.ff[form]
                )
                if result is not None:
                    return ('equivalent', result[0], form, result[2])

        # try automatic equivalences
        if 'auto_equivalence' in self.ff:
            iauto = self.ff['auto_equivalence'][i]['torsion_end_atom']
            jauto = self.ff['auto_equivalence'][j]['torsion_center_atom']
            kauto = self.ff['auto_equivalence'][k]['torsion_center_atom']
            lauto = self.ff['auto_equivalence'][l]['torsion_end_atom']
            key, flipped = self.make_canonical(
                'like_torsion', (iauto, jauto, kauto, lauto)
            )
            for form in forms:
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])

                # try wildcards, which may have numerical precidence
                # Find all the single-sided wildcards, realizing that the
                # triplet might be flipped.
                left = []
                right = []
                for key in self.ff[form]:
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
                        key, flipped = self.make_canonical(
                            'like_torsion', (left[0], jauto, kauto, lauto)
                        )
                        if key in self.ff[form]:
                            return ('automatic', key, form, self.ff[form][key])
                    else:
                        if left[0] < right[0]:
                            key, flipped = self.make_canonical(
                                'like_torsion', (left[0], jauto, kauto, lauto)
                            )
                            if key in self.ff[form]:
                                return (
                                    'automatic', key, form, self.ff[form][key]
                                )
                        else:
                            key, flipped = self.make_canonical(
                                'like_torsion',
                                (iauto, jauto, kauto, right[0])
                            )
                            if key in self.ff[form]:
                                return (
                                    'automatic', key, form, self.ff[form][key]
                                )
                elif len(right) > 0:
                    key, flipped = self.make_canonical(
                        'like_torsion', (iauto, jauto, kauto, right[0])
                    )
                    if key in self.ff[form]:
                        return ('automatic', key, form, self.ff[form][key])

                key, flipped = self.make_canonical(
                    'like_torsion', (iauto, jauto, kauto, '*')
                )
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])
                key, flipped = self.make_canonical(
                    'like_torsion', ('*', jauto, kauto, lauto)
                )
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])
                key, flipped = self.make_canonical(
                    'like_torsion', ('*', jauto, kauto, '*')
                )
                if key in self.ff[form]:
                    return ('automatic', key, form, self.ff[form][key])

        raise RuntimeError(
            'No torsion parameters for {}-{}-{}-{}'.format(i, j, k, l)
        )

    def _torsion_parameters_helper(self, i, j, k, l, section):  # noqa: E741
        """Return the torsion parameters given four atom types
        """

        # parameter directly available
        key, flipped = self.make_canonical('like_torsion', (i, j, k, l))
        if key in section:
            return (key, flipped, section[key])

        # try wildcards
        key, flipped = self.make_canonical('like_torsion', ('*', j, k, l))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical('like_torsion', (i, j, k, '*'))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical('like_torsion', ('*', j, k, '*'))
        if key in section:
            return (key, flipped, section[key])

        return None

    def oop_parameters(self, i, j, k, l, zero=False):  # noqa: E741
        """Return the oop parameters given four atoms types

        Handles equivalences and automatic equivalences and wildcards,
        with numerical precedences
        """

        forms = self.ff['terms']['out-of-plane']

        for form in forms:
            result = self._oop_parameters_helper(i, j, k, l, form)
            if result is not None:
                return ('explicit', result[0], form, result[1])

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['oop']
            jeq = self.ff['equivalence'][j]['oop']
            keq = self.ff['equivalence'][k]['oop']
            leq = self.ff['equivalence'][l]['oop']
            for form in forms:
                result = self._oop_parameters_helper(ieq, jeq, keq, leq, form)
                if result is not None:
                    return ('equivalent', result[0], form, result[1])

        # try automatic equivalences
        if 'auto_equivalence' in self.ff:
            iauto = self.ff['auto_equivalence'][i]['oop_end_atom']
            jauto = self.ff['auto_equivalence'][j]['oop_center_atom']
            kauto = self.ff['auto_equivalence'][k]['oop_end_atom']
            lauto = self.ff['auto_equivalence'][l]['oop_end_atom']
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
        key, flipped = self.make_canonical('like_oop', (i, j, k, l))
        if key in self.ff[form]:
            return (key, self.ff[form][key])

        # try wildcards
        key, flipped = self.make_canonical('like_oop', ('*', j, k, l))
        if key in self.ff[form]:
            return (key, self.ff[form][key])
        key, flipped = self.make_canonical('like_oop', (i, j, '*', l))
        if key in self.ff[form]:
            return (key, self.ff[form][key])
        key, flipped = self.make_canonical('like_oop', (i, j, k, '*'))
        if key in self.ff[form]:
            return (key, self.ff[form][key])
        key, flipped = self.make_canonical('like_oop', ('*', j, '*', l))
        if key in self.ff[form]:
            return (key, self.ff[form][key])
        key, flipped = self.make_canonical('like_oop', ('*', j, k, '*'))
        if key in self.ff[form]:
            return (key, self.ff[form][key])
        key, flipped = self.make_canonical('like_oop', (i, j, '*', '*'))
        if key in self.ff[form]:
            return (key, self.ff[form][key])
        key, flipped = self.make_canonical('like_oop', ('*', j, '*', '*'))
        if key in self.ff[form]:
            return (key, self.ff[form][key])

        return None

    def nonbond_parameters(self, i, j=None, form='nonbond(12-6)'):
        """Return the nondbond parameters given one or two atoms types i and j

        Handle equivalences
        """

        # parameter directly available
        if j is None:
            key = (i,)
        else:
            key, flipped = self.make_canonical('like_bond', (i, j))
        if key in self.ff[form]:
            return ('explicit', key, form, self.ff[form][key])

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['nonbond']
            if j is None:
                key = (ieq,)
            else:
                jeq = self.ff['equivalence'][j]['nonbond']
                key, flipped = self.make_canonical('like_bond', (ieq, jeq))
            if key in self.ff[form]:
                return ('equivalent', key, form, self.ff[form][key])

        # try automatic equivalences
        if 'auto_equivalence' in self.ff:
            iauto = self.ff['auto_equivalence'][i]['nonbond']
            if j is None:
                key = (iauto,)
            else:
                jauto = self.ff['auto_equivalence'][j]['nonbond']
                key, flipped = self.make_canonical('like_bond', (iauto, jauto))
            if key in self.ff[form]:
                return ('automatic', key, form, self.ff[form][key])

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
        result = self._angle_parameters_helper(i, j, k, self.ff['bond-bond'])
        if result is not None:
            if result[1]:
                values = {
                    'R10': b2_parameters['R0'],
                    'R20': b1_parameters['R0']
                }
            values.update(result[2])
            return ('explicit', result[0], 'bond-bond', values)

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['angle']
            jeq = self.ff['equivalence'][j]['angle']
            keq = self.ff['equivalence'][k]['angle']
            result = self._angle_parameters_helper(
                ieq, jeq, keq, self.ff['bond-bond']
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
        key, flipped = self.make_canonical('like_angle', (i, j, k))
        if key in section:
            return (key, flipped, section[key])

        # try wildcards
        key, flipped = self.make_canonical('like_angle', ('*', j, k))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical('like_angle', (i, j, '*'))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical('like_angle', ('*', j, '*'))
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
            i, j, k, l, self.ff['bond-bond_1_3']
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
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['torsion']
            jeq = self.ff['equivalence'][j]['torsion']
            keq = self.ff['equivalence'][k]['torsion']
            leq = self.ff['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.ff['bond-bond_1_3']
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
        result = self._angle_parameters_helper(i, j, k, self.ff['bond-angle'])
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
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['angle']
            jeq = self.ff['equivalence'][j]['angle']
            keq = self.ff['equivalence'][k]['angle']
            result = self._angle_parameters_helper(
                ieq, jeq, keq, self.ff['bond-angle']
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
            i, j, k, l, self.ff['angle-angle']
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
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['angle']
            jeq = self.ff['equivalence'][j]['angle']
            keq = self.ff['equivalence'][k]['angle']
            leq = self.ff['equivalence'][l]['angle']
            result = self._angle_angle_parameters_helper(
                ieq, jeq, keq, leq, self.ff['angle-angle']
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
        key, flipped = self.make_canonical('like_angle-angle', (i, j, k, l))
        if key in section:
            return (key, flipped, section[key])

        # try wildcards
        key, flipped = self.make_canonical('like_angle-angle', ('*', j, k, l))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical('like_angle-angle', (i, j, k, '*'))
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical(
            'like_angle-angle', ('*', j, k, '*')
        )
        if key in section:
            return (key, flipped, section[key])
        key, flipped = self.make_canonical(
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
            i, j, k, l, self.ff['end_bond-torsion_3']
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
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['torsion']
            jeq = self.ff['equivalence'][j]['torsion']
            keq = self.ff['equivalence'][k]['torsion']
            leq = self.ff['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.ff['end_bond-torsion_3']
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
            i, j, k, l, self.ff['middle_bond-torsion_3']
        )
        if result is not None:
            values.update(result[2])
            return ('explicit', result[0], 'middle_bond-torsion_3', values)

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['torsion']
            jeq = self.ff['equivalence'][j]['torsion']
            keq = self.ff['equivalence'][k]['torsion']
            leq = self.ff['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.ff['middle_bond-torsion_3']
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
            i, j, k, l, self.ff['angle-torsion_3']
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
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['torsion']
            jeq = self.ff['equivalence'][j]['torsion']
            keq = self.ff['equivalence'][k]['torsion']
            leq = self.ff['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.ff['angle-torsion_3']
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
            i, j, k, l, self.ff['angle-angle-torsion_1']
        )
        if result is not None:
            values.update(result[2])
            return ('explicit', result[0], 'angle-angle-torsion_1', values)

        # try equivalences
        if 'equivalence' in self.ff:
            ieq = self.ff['equivalence'][i]['torsion']
            jeq = self.ff['equivalence'][j]['torsion']
            keq = self.ff['equivalence'][k]['torsion']
            leq = self.ff['equivalence'][l]['torsion']
            result = self._torsion_parameters_helper(
                ieq, jeq, keq, leq, self.ff['angle-angle-torsion_1']
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

    def get_templates(self):
        """Return the templates dict
        """
        return self.ff['templates']


    def eex_charges(self, eex, system, configuration=None):
        """Do nothing routine since charges are handled by the increments."""
        pass

    def eex_increment(self, eex, system, configuration=None):
        """Get the charges for the structure

        If they do not exists on the structure, they are created
        using the bond increments and saved on the structure"""

        logger.debug('entering eex_increment')

        ff_name = self.name
        atoms = system['atom']
        key = f'charges_{ff_name}'
        if key in atoms:
            eex['charges'] = [*atoms[key]]
        else:
            raise RuntimeError('No charges on system!')

        logger.debug('leaving eex_increment')


    def eex_pair(self, eex, system, configuration=None):
        """Create the pair (non-bond) portion of the energy expression"""
        logger.debug('In eex_pair')
        types = self.topology['types']

        for pair_type in ('nonbond(12-6)', 'nonbond(9-6)'):
            if pair_type in self.ff['functional_forms']:
                found = True
                break
        if not found:
            raise RuntimeError('Error finding pair_type in eex_pair')

        result = eex['nonbonds'] = []
        parameters = eex['nonbond parameters'] = []
        for itype in types[1:]:
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
