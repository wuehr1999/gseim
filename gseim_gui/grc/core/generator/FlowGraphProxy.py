# Copyright 2016 Free Software Foundation, Inc.
# This file is part of GNU Radio
#
# GNU Radio Companion is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# GNU Radio Companion is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA

from __future__ import absolute_import
from six.moves import range

class FlowGraphProxy(object):  # TODO: move this in a refactored Generator

    def __init__(self, fg):
        self.orignal_flowgraph = fg

    def __getattr__(self, item):
        return getattr(self.orignal_flowgraph, item)

    def get_hier_block_io(self, direction):
        """
        Get a list of io ports for this flow graph.

        Args:
            direction: a string of 'in' or 'out'

        Returns:
            a list of dicts with: type, label, vlen, size, optional
        """
        pads = self.get_pad_sources() if direction in ('sink', 'in') else \
            self.get_pad_sinks() if direction in ('source', 'out') else []
        ports = []
        for pad in pads:
            type_param = pad.params['type']
            master = {
                'label': str(pad.params['label'].get_evaluated()),
                'type': str(pad.params['type'].get_evaluated()),
                'vlen': str(pad.params['vlen'].get_value()),
                'size':  type_param.options.attributes[type_param.get_value()]['size'],
                'cpp_size':  type_param.options.attributes[type_param.get_value()]['cpp_size'],
                'optional': bool(pad.params['optional'].get_evaluated()),
            }
            num_ports = pad.params['num_streams'].get_evaluated()
            if num_ports > 1:
                for i in range(num_ports):
                    clone = master.copy()
                    clone['label'] += str(i)
                    ports.append(clone)
            else:
                ports.append(master)
        return ports

    def get_pad_sources(self):
        """
        Get a list of pad source blocks sorted by id order.

        Returns:
            a list of pad source blocks in this flow graph
        """
        pads = [b for b in self.get_enabled_blocks() if b.key == 'pad_source']
        return sorted(pads, key=lambda x: x.name)

    def get_pad_sinks(self):
        """
        Get a list of pad sink blocks sorted by id order.

        Returns:
            a list of pad sink blocks in this flow graph
        """
        pads = [b for b in self.get_enabled_blocks() if b.key == 'pad_sink']
        return sorted(pads, key=lambda x: x.name)

def get_hier_block_io(flow_graph, direction, domain=None):
    """
    Get a list of io ports for this flow graph.

    Returns a list of dicts with: type, label, vlen, size, optional
    """
    pads = flow_graph.get_pad_sources() if direction in ('sink', 'in') else \
        flow_graph.get_pad_sinks() if direction in ('source', 'out') else []
    ports = []
    for pad in pads:
        type_param = pad.params['type']
        master = {
            'label': str(pad.params['label'].get_evaluated()),
            'type': str(pad.params['type'].get_evaluated()),
            'vlen': str(pad.params['vlen'].get_value()),
            'size':  type_param.options.attributes[type_param.get_value()]['size'],
            'optional': bool(pad.params['optional'].get_evaluated()),
        }
        num_ports = pad.params['num_streams'].get_evaluated()
        if num_ports > 1:
            for i in range(num_ports):
                clone = master.copy()
                clone['label'] += str(i)
                ports.append(clone)
        else:
            ports.append(master)
    if domain is not None:
        ports = [p for p in ports if p.domain == domain]
    return ports
