# -*- coding: utf-8 -*-

"""A node in a flowchart


"""

import abc
import json
import logging
import seamm
import seamm_util  # MUST come after seamm
from seamm_util.printing import FormattedText as __
import seamm_util.printing as printing
import numpy as np
import os.path
import pandas
import uuid

logger = logging.getLogger(__name__)
job = printing.getPrinter()


class Node(abc.ABC):

    def __init__(self, flowchart=None, title='', extension=None):
        """Initialize a node

        Keyword arguments:
        """

        self._uuid = uuid.uuid4().int
        self.parent = None
        self.flowchart = flowchart
        self._title = title
        self._description = ''
        self._id = None
        self.extension = extension
        self._visited = False

        self.parameters = None  # Object containing control parameters

        self.x = None
        self.y = None
        self.w = None
        self.h = None

        # Set up our formatter for printing
        self.formatter = logging.Formatter(fmt='{message:s}', style='{')

    def __hash__(self):
        """Make iterable!"""
        return self._uuid

    @property
    def uuid(self):
        """The uuid of the node"""
        return self._uuid

    @property
    def title(self):
        """The title to display"""
        return self._title

    @title.setter
    def title(self, value):
        self._title = value

    @property
    def tag(self):
        """The string representation of the uuid of the node"""
        return 'node=' + str(self._uuid)

    @property
    def directory(self):
        """Return the directory we should write output to"""
        return os.path.join(self.flowchart.root_directory, *self._id)

    @property
    def visited(self):
        """Whether this node has been visited in a traversal"""
        return self._visited

    @visited.setter
    def visited(self, value):
        self._visited = bool(value)

    @property
    def description(self):
        """A textual description of this node"""
        return self._description

    @description.setter
    def description(self, value):
        self._description = value

    @property
    def indent(self):
        length = len(self._id)
        if length <= 1:
            return ''
        if length > 2:
            result = (length - 2) * '  .' + '   '
        else:
            result = '   '
        return result

    @property
    def header(self):
        """A printable header for this section of output"""
        return (
            'Step {}: {}  {}'.format(
                '.'.join(str(e) for e in self._id), self.title, self.version
            )
        )

    def set_uuid(self):
        self._uuid = uuid.uuid4().int

        # Need to correct all edges to other nodes
        raise NotImplementedError('set_uuid not implemented yet!')

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            return self.next()

    def reset_id(self):
        """Reset the id for node"""
        self._id = None

    def get_gui_data(self, key, gui=None):
        """Return an element from the GUI dictionary"""
        if gui is None:
            return self._gui_data[seamm.Flowchart.graphics][key]
        else:
            return self._gui_data[gui][key]

    def set_gui_data(self, key, value, gui=None):
        """Set an element of the GUI dictionary"""
        if gui is None:
            if seamm.Flowchart.graphics not in self._gui_data:
                self._gui_data[seamm.Flowchart.graphics] = {}
            self._gui_data[seamm.Flowchart.graphics][key] = value
        else:
            if gui not in self._gui_data:
                self._gui_data[gui] = {}
            self._gui_data[gui][key] = value

    def connections(self):
        """Return a list of all the incoming and outgoing edges
        for this node, giving the anchor points and other node
        """

        result = self.flowchart.edges(self)
        return result

    def remove_edge(self, edge):
        """Remove a given edge, or all edges if 'all' is given
        """

        if isinstance(edge, str) and edge == 'all':
            for direction, obj in self.connections():
                self.remove_edge(obj)
        else:
            self.flowchart.graph.remove_edge(
                edge.node1, edge.node2, edge.edge_type, edge.edge_subtype
            )

    def description_text(self, P=None):
        """Prepare information about what this node will do
        """
        return (
            'This node has no specific description. '
            "Override the method 'description_text' "
            'to provide the description.'
        )

    def describe(self):
        """Write out information about what this node will do
        """

        self.visited = True

        # The 'step' line
        job.job('')
        job.job(__(self.header, indent=self.indent))

        # and rest of the description
        if self.parameters:
            P = self.parameters.values_to_dict()
            text = self.description_text(P)
            job.job(__(text, **P, indent=self.indent + '    '))

        next_node = self.next()

        if next_node is None or next_node.visited:
            return None
        else:
            return next_node

    def run(self, printer=None):
        """Do whatever we need to do! The base class does nothing except
        return the next node.
        """

        # Create the directory for writing output and files
        os.makedirs(self.directory, exist_ok=True)

        if printer is not None:
            # Setup up the printing for this step
            self.setup_printing(printer)

            printer.important(self.header)

        next_node = self.next()
        if next_node:
            logger.debug('returning next_node: {}'.format(next_node))
        else:
            logger.debug('returning next_node: None')

        return next_node

    def next(self):
        """Return the next node in the flow
        """

        for edge in self.flowchart.edges(self, direction='out'):
            if edge.edge_subtype == 'next':
                logger.debug('Next node is: {}'.format(edge.node2))
                return edge.node2

        logger.debug('Reached the end of the flowchart')
        return None

    def previous(self):
        """Return the previous node in the flow
        """

        for edge in self.flowchart.edges(self, direction='in'):
            if edge.edge_type == 'execution' and \
               edge.edge_subtype == 'next':
                return edge.node1

        return None

    def get_input(self):
        """Return the input from this subnode, usually used for
        building up the input for the executable."""

        return ''

    def to_json(self):
        return json.dumps(self.to_dict(), cls=seamm_util.JSONEncoder)

    def to_dict(self):
        """serialize this object and everything it contains as a dict"""
        data = {
            'item': 'object',
            'module': self.__module__,
            'class': self.__class__.__name__,
            'extension': self.extension
        }
        data['attributes'] = {}
        for key in self.__dict__:
            if key == 'flowchart':
                continue
            if key == 'parent':
                continue
            if key == 'formatter':
                continue
            if 'flowchart' in key:
                # Have a subflowchart!
                data[key] = self.__dict__[key].to_dict()
            else:
                data['attributes'][key] = self.__dict__[key]
        return data

    def from_dict(self, data):
        """un-serialize object and everything it contains from a dict"""
        if data['item'] != 'object':
            raise RuntimeError('The data for restoring the object is invalid')
        if data['class'] != self.__class__.__name__:
            raise RuntimeError(
                'Trying to restore a {}'.format(self.__class__.__name__) +
                ' from data for a {}'.format(data['class'])
            )
        for key in data:
            if key == 'attributes':
                attributes = data['attributes']
                for key in attributes:
                    self.__dict__[key] = attributes[key]
            elif 'flowchart' in key:
                self.__dict__[key].from_dict(data[key])

    def default_edge_subtype(self):
        """Return the default subtype of the edge. Usually this is 'next'
        but for nodes with two or more edges leaving them, such as a loop, this
        method will return an appropriate default for the current edge. For
        example, by default the first edge emanating from a loop-node is the
        'loop' edge; the second, the 'exit' edge.

        A return value of 'too many' indicates that the node exceeds the number
        of allowed exit edges.
        """

        # how many outgoing edges are there?
        n_edges = len(self.flowchart.edges(self, direction='out'))

        logger.debug('node.default_edge_subtype, n_edges = {}'.format(n_edges))

        if n_edges == 0:
            return ""
        else:
            return "too many"

    def analyze(self, indent='', **kwargs):
        """Analyze the output of the calculation
        """
        return

    def get_value(self, variable_or_value):
        """Return the value of the workspace variable is <variable_or_value>
        is the name of a variable. Otherwise, simply return the value of
        <variable_or_value>.

        This provides a convenient way to handle both values and variables
        in widgets. A reference to a variable is $<name> or ${name}, and is
        replaced by the contents of the variable. If the text is not a
        reference to a variable then the value passed in is returned
        unchanged.
        """

        return seamm.flowchart_variables.value(variable_or_value)

    def get_variable(self, variable):
        """Get the value of a variable, which must exist
        """

        return seamm.flowchart_variables.get_variable(variable)

    def set_variable(self, variable, value):
        """Set the value of a variable in the workspace. The name of the
        variable maybe a plain string, or be $<name> or ${<name>}
        """

        seamm.flowchart_variables.set_variable(variable, value)

    def variable_exists(self, variable):
        """Return whether a varable exists in the workspace
        """

        return seamm.flowchart_variables.exists(variable)

    def delete_variable(self, variable):
        """Delete a variable in the workspace
        """

        seamm.flowchart_variables.delete(variable)

    def setup_printing(self, printer):
        """Establish the handlers for printing as controlled by
        options
        """

        # First remove an existing handlers
        self.close_printing(printer)

        # A handler for stdout
        console_handler = logging.StreamHandler()
        console_handler.setLevel(printing.JOB)
        console_handler.setFormatter(self.formatter)
        printer.addHandler(console_handler)

        # A handler for the file
        file_handler = logging.FileHandler(
            os.path.join(self.directory, 'step.out'), delay=True
        )
        file_handler.setLevel(printing.NORMAL)
        file_handler.setFormatter(self.formatter)
        printer.addHandler(file_handler)

    def close_printing(self, printer):
        """Close the handlers for printing, so that buffers are
        flushed, files closed, etc.
        """
        for handler in printer.handlers:
            printer.removeHandler(handler)

    def job_output(self, text):
        """Temporary!"""
        job.job(text)

    def store_results(
        self, data={}, properties=None, results=None, create_tables=True
    ):
        """Store results in variables and tables, as requested

        Keywords:

        properties (dict): a dictionary of properties
        results (dict): a dictionary of results from the calculation
        create_tables (bool): whether to create tables as needed

        Each item in 'results' is itself a dictionary. If the following keys
        are in the dictionary, the appropriate action is taken:

        'variable' -- is the name of a variable to store the result in
        'table' -- the name of the table, and
        'column' -- is the column name for the result in the table.
        """

        for key, value in results.items():
            # Check for storing in a variable
            if 'variable' in value:
                self.set_variable(value['variable'], data[key])

            # and table
            if 'table' in value:
                tablename = value['table']
                column = value['column']
                # Does the table exist?
                if not self.variable_exists(tablename):
                    if create_tables:
                        table = pandas.DataFrame()
                        self.set_variable(
                            tablename, {
                                'type': 'pandas',
                                'table': table,
                                'defaults': {},
                                'loop index': False,
                                'current index': 0
                            }
                        )
                    else:
                        raise RuntimeError(
                            "Table '{}' does not exist.".format(tablename)
                        )

                table_handle = self.get_variable(tablename)
                table = table_handle['table']

                # create the column as needed
                if column not in table.columns:
                    kind = properties[key]['type']
                    if kind == 'boolean':
                        default = False
                    elif kind == 'integer':
                        default = 0
                    elif kind == 'float':
                        default = np.nan
                    else:
                        default = ''

                    table_handle['defaults'][column] = default
                    table[column] = default

                # and put the value in (finally!!!)
                row_index = table_handle['current index']
                table.at[row_index, column] = data[key]
