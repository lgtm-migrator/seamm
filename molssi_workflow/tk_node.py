# -*- coding: utf-8 -*-

import molssi_workflow
import tkinter as tk
"""A graphical node using Tk on a canvas"""

anchor_points = {
    's': (0, 1),
    'sse': (0.25, 1),
    'se': (0.5, 1),
    'ese': (0.5, 0.75),
    'e': (0.5, 0.5),
    'ene': (0.5, 0.25),
    'ne': (0.5, 0),
    'nne': (0.25, 0),
    'n': (0, 0),
    'nnw': (-0.25, 0),
    'nw': (-0.5, 0),
    'wnw': (-0.5, 0.25),
    'w': (-0.5, 0.5),
    'wsw': (-0.5, 0.75),
    'sw': (-0.5, 1),
    'ssw': (-0.25, 1)
}


class TkNode(object):
    """The node_class is the class of the 'real' node that this
    class is the Tk graphics partner for
    """

    node_class = molssi_workflow.Node

    def __init__(self, node=None, canvas=None, x=None, y=None, w=None, h=None):
        """Initialize a node

        Keyword arguments:
        """

        self._node = node
        self.toplevel = None
        self.canvas = canvas
        self._tmp_gui_data = {
            "x": None,
            "y": None,
            "w": None,
            "h": None,
        }
        if x is not None:
            self.x = x
        if y is not None:
            self.y = y
        if w is not None:
            self.w = w
        if h is not None:
            self.h = h

        self._border = None
        self.title_label = None
        self._selected = False
        self.popup_menu = None
        self._tmp = None
        self.previous_grab = None

    def __hash__(self):
        """Make iterable!"""
        return self.node.uuid

    @property
    def uuid(self):
        """The uuid of the node"""
        return self.node.uuid

    @property
    def title(self):
        """The title to display"""
        return self.node.title

    @title.setter
    def title(self, value):
        self.node.title = value
        if self.title_label is not None:
            self.canvas.itemconfigure(self.title_label, text=value)

    @property
    def tag(self):
        """The string representation of the uuid of the node"""
        return self.node.tag

    @property
    def workflow(self):
        """The workflow object"""
        return self.node.workflow

    @workflow.setter
    def workflow(self, value):
        """The workflow object"""
        self.node.workflow = value

    @property
    def x(self):
        """The x coordinate of the top center of the node"""
        if self._node is None:
            return self._tmp_gui_data['x']
        else:
            return self._node.get_gui_data('x', gui='Tk')

    @x.setter
    def x(self, value):
        if self._node is None:
            self._tmp_gui_data['x'] = value
        else:
            self._node.set_gui_data('x', value, gui='Tk')

    @property
    def y(self):
        """The y coordinate of the top center of the node"""
        if self._node is None:
            return self._tmp_gui_data['y']
        else:
            return self._node.get_gui_data('y', gui='Tk')

    @y.setter
    def y(self, value):
        if self._node is None:
            self._tmp_gui_data['y'] = value
        else:
            self._node.set_gui_data('y', value, gui='Tk')

    @property
    def w(self):
        """The width of the graphical node"""
        if self._node is None:
            return self._tmp_gui_data['w']
        else:
            return self._node.get_gui_data('w', gui='Tk')

    @w.setter
    def w(self, value):
        if self._node is None:
            self._tmp_gui_data['w'] = value
        else:
            self._node.set_gui_data('w', value, gui='Tk')

    @property
    def h(self):
        """The height of the graphical node"""
        if self._node is None:
            return self._tmp_gui_data['h']
        else:
            return self._node.get_gui_data('h', gui='Tk')

    @h.setter
    def h(self, value):
        if self._node is None:
            self._tmp_gui_data['h'] = value
        else:
            self._node.set_gui_data('h', value, gui='Tk')

    @property
    def node(self):
        """The non-graphical node we represent"""
        return self._node

    @node.setter
    def node(self, node):
        if self._node is None:
            self._node = node
            for key in self._tmp_gui_data:
                node.set_gui_data(key, self._tmp_gui_data[key],
                                  gui='Tk')
            self._tmp_gui_data = {}
        else:
            self._node = node

    def set_uuid(self):
        self.node.set_uuid()

    def connections(self):
        """Return a list of all the incoming and outgoing edges
        for this node, giving the anchor points and other node
        """

        return self.node.connections()

    @property
    def selected(self):
        """Whether I am selected or not"""
        return self._selected

    @selected.setter
    def selected(self, value):
        self._selected = value
        if value:
            self.canvas.itemconfigure(self.border, outline='red')
        else:
            self.canvas.itemconfigure(self.border, outline='black')

    @property
    def canvas(self):
        """The canvas for drawing the node"""
        return self._canvas

    @canvas.setter
    def canvas(self, value):
        if value:
            self.toplevel = value.winfo_toplevel()
        self._canvas = value

    @property
    def border(self):
        """The border of the picture in the flowchart"""
        return self._border

    @border.setter
    def border(self, value):
        self._border = value

    def anchor_point(self, anchor="all"):
        """Where the anchor points are located. If "all" is given
        a dictionary of all points is returned"""

        if anchor == "all":
            result = []
            for pt in anchor_points:
                a, b = anchor_points[pt]
                result.append((pt, int(self.x + a * self.w),
                               int(self.y + b * self.h)))
            return result

        if anchor in anchor_points:
            a, b = anchor_points[anchor]
            return (int(self.x + a * self.w), int(self.y + b * self.h))

        raise NotImplementedError(
            "anchor position '{}' not implemented".format(anchor))

    def check_anchor_points(self, x, y, halo):
        """If the position x, y is within halo or one of the anchor points
        activate the point and return the name of the anchor point
        """

        points = []
        for direction, edge in self.connections():
            if direction == 'out':
                points.append(edge.gui_object['start_point'])
            else:
                points.append(edge.gui_object['end_point'])

        for point, x0, y0 in self.anchor_point():
            if x >= x0 - halo and x <= x0 + halo and \
               y >= y0 - halo and y <= y0 + halo:
                if point in points:
                    return None
                else:
                    return point
        return None

    def draw(self):
        """Draw the node on the given canvas, making it visible"""

        # the outline
        x0 = self.x - self.w / 2
        x1 = x0 + self.w
        y0 = self.y
        y1 = y0 + self.h
        self._border = self.canvas.create_rectangle(
            x0,
            y0,
            x1,
            y1,
            tags=[self.tag, 'type=outline'],
            fill='white',
        )

        # the label in the middle
        x0 = self.x
        y0 = self.y + self.h / 2
        self.title_label = self.canvas.create_text(
            x0, y0, text=self.title, tags=[self.tag, 'type=title'])

    def undraw(self):
        """Remove all of our visual components
        """

        self.canvas.delete(self.tag)

    def move(self, deltax, deltay):
        if self._tmp is None:
            self._tmp = self.connections()

        self.x += deltax
        self.y += deltay

        self.canvas.move(self.tag, deltax, deltay)

        for connection in self._tmp:
            direction, edge = connection
            edge.gui_object.move()

    def end_move(self, deltax, deltay):
        self.move(deltax, deltay)
        self._x0 = None
        self._y0 = None
        self._tmp = None

    def right_click(self, event):
        """Do whatever needs to be done for a right-click on this
        item in the flowchart.

        Subclasses should override this as appropriate! The menu
        created in this base method is accessible in subclasses
        which should make it relatively easy to override.
        """

        if self.popup_menu is not None:
            self.popup_menu.destroy()

        self.popup_menu = tk.Menu(self.canvas, tearoff=0)
        self.popup_menu.add_command(
            label="Delete", command=lambda: self.workflow.remove_node(self))

        if type(self) is molssi_workflow.tk_node.TkNode:
            self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def double_click(self, event):
        """Do whatever needs to be done for a double-click on this
        item in the flowchart.

        Subclasses should override this as appropriate!
        """

        self.edit()

    def activate(self):
        """Add active handles at the anchor points and change the
        cursor
        """

        self.canvas.delete(self.tag + ' && type=anchor')
        for pt, x, y in self.anchor_point("all"):
            x0 = x - 2
            y0 = y - 2
            x1 = x + 2
            y1 = y + 2
            self.canvas.create_oval(
                x0,
                y0,
                x1,
                y1,
                fill='red',
                outline='red',
                tags=[self.tag, 'type=anchor', 'anchor=' + pt])

    def deactivate(self):
        """Remove the decorations indicate the anchor points
        """

        self.canvas.delete(self.tag + ' && type=anchor')
        self.canvas.delete(self.tag + ' && type=active_anchor')

    def is_inside(self, x, y, halo=0):
        """Return a boolean indicating whether the point x, y is inside
        this node, using halo as a size around the point
        """
        if x < self.x - self.w / 2 - halo:
            return False
        if x > self.x + self.w / 2 + halo:
            return False

        if y < self.y - halo:
            return False
        if y > self.y + self.h + halo:
            return False

        return True

    def activate_anchor_point(self, point, halo):
        """Put a marker on the anchor point to indicate it is
        active
        """

        x, y = self.anchor_point(point)
        self.canvas.create_oval(
            x - halo,
            y - halo,
            x + halo,
            y + halo,
            fill='red',
            outline='red',
            tags=[self.tag, 'type=active_anchor', 'anchor=' + point])

    def remove_edge(self, edge):
        """Remove a given edge, or all edges if 'all' is given
        """

        if isinstance(edge, str) and edge == 'all':
            for direction, obj in self.connections():
                self.remove_edge(obj)
        else:
            self.canvas.delete(edge.gui_object.tag())
            self.node.remove_edge(edge)

    def edit(self):
        """Do-nothing base class method"""
        pass

    def to_dict(self):
        """Serialize to a dict"""
        data = {
            'x': self._x,
            'y': self._y,
            'w': self._w,
            'h': self._h,
        }

        return data
