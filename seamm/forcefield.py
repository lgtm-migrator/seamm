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
        self.ff['terms'] = {}
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
                    self.ff['terms'][term].append(fform)
                else:
                    self.ff['terms'][term] = [fform]

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
