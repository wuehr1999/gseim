"""
Copyright 2008-2015 Free Software Foundation, Inc.
This file is part of GNU Radio

GNU Radio Companion is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

GNU Radio Companion is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
"""

from __future__ import absolute_import

import collections
import itertools
import copy

import six
from six.moves import range
import re

import ast

from ._flags import Flags

from ..base import Element
from ..utils.descriptors import lazy_property

def _get_elem(iterable, key):
    items = list(iterable)
    for item in items:
        if item.key == key:
            return item
    return ValueError('Key "{}" not found in {}.'.format(key, items))

class Block(Element):

    is_block = True

    key = ''
    label = ''
    category = ''
    flags = Flags('')
    documentation = {'': ''}

    value = None
    asserts = []

    parameters_data = []
    inputs_data = []
    outputs_data = []

    extra_data = {}
    loaded_from = '(unknown)'

    def __init__(self, parent):
        """Make a new block from nested data."""
        super(Block, self).__init__(parent)
        param_factory = self.parent_platform.make_param
        port_factory = self.parent_platform.make_port

        self.params = collections.OrderedDict(
            (data['id'], param_factory(parent=self, **data))
            for data in self.parameters_data
        )
        if self.key == 'options':
            self.params['id'].hide = 'part'

        self.sinks = [port_factory(parent=self, **params) for params in self.inputs_data]
        self.sources = [port_factory(parent=self, **params) for params in self.outputs_data]

        self.active_sources = []  # on rewrite
        self.active_sinks = []  # on rewrite

        self.states = {'state': True, 'bus_source': False, 'bus_sink': False, 'bus_structure': None}

        self.current_bus_structure = {'source': None, 'sink': None}

    def get_bus_structure(self, direction):
        if direction == 'source':
            bus_structure = self.bus_structure_source
        else:
            bus_structure = self.bus_structure_sink

        if not bus_structure:
            return None

        try:
            clean_bus_structure = self.evaluate(bus_structure)
            return clean_bus_structure
        except:
            return None

    # region Rewrite_and_Validation
    def rewrite(self):
        """
        Add and remove ports to adjust for the nports.
        """
        Element.rewrite(self)

        def rekey(ports):
            """Renumber non-message/message ports"""
            domain_specific_port_index = collections.defaultdict(int)
            for port in ports:
                if not port.key.isdigit():
                    continue
                domain = port.domain
                port.key = str(domain_specific_port_index[domain])
                domain_specific_port_index[domain] += 1

        # Adjust nports
        for ports in (self.sources, self.sinks):
            self._rewrite_nports(ports)
            rekey(ports)

        self.update_bus_logic()
        # disconnect hidden ports
        self.parent_flowgraph.disconnect(*[p for p in self.ports() if p.hidden])

        self.active_sources = [p for p in self.sources if not p.hidden]
        self.active_sinks = [p for p in self.sinks if not p.hidden]

    def update_bus_logic(self):
        
        for direc in {'source','sink'}:
            if direc == 'source':
                ports = self.sources
                ports_gui = self.filter_bus_port(self.sources)
                bus_state = self.bus_source
            else:
                ports = self.sinks
                ports_gui = self.filter_bus_port(self.sinks)  
                bus_state = self.bus_sink

            # Remove the bus ports
            removed_bus_ports = []
            removed_bus_connections = []
            if 'bus' in map(lambda a: a.dtype, ports):
                for port in ports_gui:
                    for c in self.parent_flowgraph.connections:
                        if port is c.source_port or port is c.sink_port:
                            removed_bus_ports.append(port)
                            removed_bus_connections.append(c)
                    ports.remove(port)

            if (bus_state):
                struct = self.form_bus_structure(direc)
                self.current_bus_structure[direc] = struct

                # Hide ports that are not part of the bus structure
                #TODO: Blocks where it is desired to only have a subset 
                # of ports included in the bus still has some issues
                for idx, port in enumerate(ports):
                    if any([idx in bus for bus in self.current_bus_structure[direc]]):
                        if (port.stored_hidden_state is None):
                            port.stored_hidden_state = port.hidden
                            port.hidden = True

                # Add the Bus Ports to the list of ports
                for i in range(len(struct)):
                    # self.sinks = [port_factory(parent=self, **params) for params in self.inputs_data]
                    port = self.parent.parent.make_port(self,direction=direc,id=str(len(ports)),label='bus',dtype='bus',bus_struct=struct[i])
                    ports.append(port)

                    for (saved_port, connection) in zip(removed_bus_ports, removed_bus_connections):
                        if port.key == saved_port.key:
                            self.parent_flowgraph.connections.remove(connection)
                            if saved_port.is_source:
                                connection.source_port = port 
                            if saved_port.is_sink:
                                connection.sink_port = port 
                            self.parent_flowgraph.connections.add(connection)

            else:
                self.current_bus_structure[direc] = None

                # Re-enable the hidden property of the ports
                for port in ports:
                    port.hidden = port.stored_hidden_state
                    port.stored_hidden_state = None

    def _rewrite_nports(self, ports):
        for port in ports:
            if hasattr(port, 'master_port'):  # Not a master port and no left-over clones
                port.dtype = port.master_port.dtype
                continue
            nports = port.multiplicity
            for clone in port.clones[nports-1:]:
                # Remove excess connections
                self.parent_flowgraph.disconnect(clone)
                port.remove_clone(clone)
                ports.remove(clone)
            # Add more cloned ports
            for j in range(1 + len(port.clones), nports):
                clone = port.add_clone()
                ports.insert(ports.index(port) + j, clone)

    def validate(self):
        """
        Validate this block.
        Call the base class validate.
        Evaluate the checks: each check must evaluate to True.
        """
        Element.validate(self)
        self._run_asserts()
        self._validate_var_value()

    def _run_asserts(self):
        """Evaluate the checks"""
        for expr in self.asserts:
            try:
                if not self.evaluate(expr):
                    self.add_error_message('Assertion "{}" failed.'.format(expr))
            except Exception:
                self.add_error_message('Assertion "{}" did not evaluate.'.format(expr))

    def _validate_var_value(self):
        """or variables check the value (only if var_value is used)"""
        if self.is_variable and self.value != 'value':
            try:
                self.parent_flowgraph.evaluate(self.value, local_namespace=self.namespace)
            except Exception as err:
                self.add_error_message('Value "{}" cannot be evaluated:\n{}'.format(self.value, err))
    # endregion

    # region Properties

    def __str__(self):
        return 'Block - {} - {}({})'.format(self.name, self.label, self.key)

    def __repr__(self):
        try:
            name = self.name
        except Exception:
            name = self.key
        return 'block[' + name + ']'

    @property
    def name(self):
        return self.params['id'].value

    @lazy_property
    def is_variable(self):
        return bool(self.value)

    @lazy_property
    def is_import(self):
        return self.key == 'import'

    @lazy_property
    def is_snippet(self):
        return self.key == 'snippet'

    @property
    def comment(self):
        return self.params['comment'].value

    @property
    def state(self):
        return 'enabled'

    @state.setter
    def state(self, value):
        """Sets the state for the block."""
        self.states['state'] = value

    # Enable/Disable Aliases
    @property
    def enabled(self):
        """Get the enabled state of the block"""
        return self.state != 'disabled'

    @property
    def bus_sink(self):
        """Gets the block's current Toggle Bus Sink state."""
        return self.states['bus_sink']

    @bus_sink.setter
    def bus_sink(self, value):
        """Sets the Toggle Bus Sink state for the block."""
        self.states['bus_sink'] = value

    @property
    def bus_source(self):
        """Gets the block's current Toggle Bus Sink state."""
        return self.states['bus_source']

    @bus_source.setter
    def bus_source(self, value):
        """Sets the Toggle Bus Source state for the block."""
        self.states['bus_source'] = value

    @property
    def bus_structure_source(self):
        """Gets the block's current source bus structure."""
        try:  
            bus_structure = self.params['bus_structure_source'].value or None
        except:
            bus_structure = None
        return bus_structure

    @property
    def bus_structure_sink(self):
        """Gets the block's current source bus structure."""
        try:  
            bus_structure = self.params['bus_structure_sink'].value or None
        except:
            bus_structure = None
        return bus_structure

    # endregion

    def ports(self):
        return itertools.chain(self.sources, self.sinks)

    def active_ports(self):
        return itertools.chain(self.active_sources, self.active_sinks)

    def children(self):
        return itertools.chain(six.itervalues(self.params), self.ports())

    # Access

    def get_sink(self, key):
        return _get_elem(self.sinks, key)

    def get_source(self, key):
        return _get_elem(self.sources, key)

    # Resolve
    @property
    def namespace(self):
        return {key: param.get_evaluated() for key, param in six.iteritems(self.params)}

    def evaluate(self, expr):
        return self.parent_flowgraph.evaluate(expr, self.namespace)

    # Import/Export Methods
    def export_data(self):
        """
        Export this block's params to nested data.

        Returns:
            a nested data odict
        """
        data = collections.OrderedDict()
        if self.key != 'options':
            data['name'] = self.name
            data['id'] = self.key
        data['parameters'] = collections.OrderedDict(sorted(
            (param_id, param.value) for param_id, param in self.params.items()
            if (param_id != 'id' or self.key == 'options')
        ))
        data['states'] = collections.OrderedDict(sorted(self.states.items()))

        return data

    def import_data(self, name, states, parameters, **_):
        """
        Import this block's params from nested data.
        Any param keys that do not exist will be ignored.
        Since params can be dynamically created based another param,
        call rewrite, and repeat the load until the params stick.
        """
        self.params['id'].value = name
        self.states.update(states)

        def get_hash():
            return hash(tuple(hash(v) for v in self.params.values()))

        pre_rewrite_hash = -1
        while pre_rewrite_hash != get_hash():
            for key, value in six.iteritems(parameters):
                try:
                    self.params[key].set_value(value)
                except KeyError:
                    continue
            # Store hash and call rewrite
            pre_rewrite_hash = get_hash()
            self.rewrite()
            
    # Controller Modify
    def filter_bus_port(self, ports):
        buslist = [p for p in ports if p.dtype == 'bus']
        return buslist or ports

    def type_controller_modify(self, direction):
        """
        Change the type controller.

        Args:
            direction: +1 or -1

        Returns:
            true for change
        """
        changed = False
        type_param = None
        for param in filter(lambda p: p.is_enum(), self.get_params()):
            children = self.get_ports() + self.get_params()
            # Priority to the type controller
            if param.get_key() in ' '.join(map(lambda p: p._type, children)): type_param = param
            # Use param if type param is unset
            if not type_param:
                type_param = param
        if type_param:
            # Try to increment the enum by direction
            try:
                keys = type_param.get_option_keys()
                old_index = keys.index(type_param.get_value())
                new_index = (old_index + direction + len(keys)) % len(keys)
                type_param.set_value(keys[new_index])
                changed = True
            except:
                pass
        return changed

    def form_bus_structure(self, direc):
        if direc == 'source':
            ports = self.sources
            bus_structure = self.get_bus_structure('source')
        else:
            ports = self.sinks
            bus_structure = self.get_bus_structure('sink')

        struct = [range(len(ports))]
        # struct = list(range(len(ports)))
        #TODO for more complicated port structures, this code is needed but not working yet
        if any([p.multiplicity for p in ports]):
            structlet = []
            last = 0
            # group the ports with > n inputs together on the bus
            cnt = 0
            idx = 0
            for p in ports:
                if cnt > 0:
                    cnt -= 1
                    continue

                if p.multiplicity > 1:
                    cnt = p.multiplicity-1
                    structlet.append([idx+j for j in range(p.multiplicity)])
                else:
                    structlet.append([idx])

            struct = structlet
        if bus_structure:
            struct = bus_structure

        self.current_bus_structure[direc] = struct
        return struct

    def bussify(self, direc):
        if direc == 'source':
            ports = self.sources
            ports_gui = self.filter_bus_port(self.sources)
            self.bus_structure = self.get_bus_structure('source')
            self.bus_source = not self.bus_source
        else:
            ports = self.sinks
            ports_gui = self.filter_bus_port(self.sinks)
            self.bus_structure = self.get_bus_structure('sink')
            self.bus_sink = not self.bus_sink

        # Disconnect all the connections when toggling the bus state
        for port in ports:
            l_connections = list(port.connections())
            for connect in l_connections:
                self.parent.remove_element(connect)

        self.update_bus_logic()
