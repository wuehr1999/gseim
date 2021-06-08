"""
Copyright (C) 2021 - Mahesh Patil <mbpatil@ee.iitb.ac.in>
This file is part of GSEIM.

GSEIM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import yaml
import sys
import os
import gutils_gseim as gu

def treat_virtual(d):

    for block in d['blocks']:
        if 'name' in block.keys():
            old_name = block['name']
            if old_name.startswith('dummy'):
                if 'parameters' in block.keys():
                    node_name = block['parameters']['n']
                    if node_name == 'none':
                        print('treat_virtual: node_name: none not allowed. Halting...',
                          flush=True)
                        sys.exit()
                else:
                    print('treat_virtual: parameters not found. Halting...',
                      flush=True)
                    sys.exit()

                l_1 = old_name.split('$')
                new_name = l_1[0] + '$' + node_name
                block['name'] = new_name

                for c in d['connections']:
                    if c[0] == old_name: c[0] = new_name
                    if c[2] == old_name: c[2] = new_name

                for k, v in d['outvars'].items():
                    if v[0] == 'connection':
                        if v[1][0] == old_name:
                            d['outvars'][k] = ['connection', [new_name, v[1][1], v[1][2], v[1][3]]]
                        if v[1][2] == old_name:
                            d['outvars'][k] = ['connection', [v[1][0], v[1][1], new_name, v[1][3]]]

    check_dummy_1(d)
    return d

def check_dummy_1(d):

#   Check that, for each dummy_source, there is a dummy_sink (and vice versa).

    s_blocks = set()

    for block in d['blocks']:
        if 'name' in block.keys():
            if block['name'].startswith('dummy'): s_blocks.add(block['name'])

    for s1 in s_blocks:
        if s1.startswith('dummy_source'):
            s2 = s1.replace('source', 'sink')
            if s2 not in s_blocks:
                print('check_dummy_1: complement of',
                  s1, ' not found. Halting...', flush=True)
                sys.exit()

        if s1.startswith('dummy_sink'):
            s2 = s1.replace('sink', 'source')
            if s2 not in s_blocks:
                print('check_dummy_1: complement of',
                  s1, ' not found. Halting...', flush=True)
                sys.exit()

#   check that there is only one dummy_sink for any node name:

    l_dummy_sinks = []

    for block in d['blocks']:
        if 'name' in block.keys():
            if block['name'].startswith('dummy_sink'):
                l_dummy_sinks.append(block['name'])

    for k in l_dummy_sinks:
        if l_dummy_sinks.count(k) != 1:
            print('check_dummy_1:', k, 'appears more than once.',
              'Halting...', flush=True)
            sys.exit()

def get_unique_names_1(l_in):
#   Note that l_in is a list of lists (GRC connections)
#   example:
#    [
#     ['gs_add$0', '0', 'gs_add$1', '0'],
#     ['gs_add$1', '0', 'pad_sink$0', '0']
#    ]

    l_out = []

    for l0 in l_in:
        s1 = l0[0]
        if s1 not in l_out:
            l_out.append(s1)
        s1 = l0[2]
        if s1 not in l_out:
            l_out.append(s1)
    return l_out

class lib_entity:
    def __init__(self, name, filename):
        self.name = name
        self.pad_in_names = []
        self.pad_out_names = []
        self.gseim_in_names = []
        self.gseim_out_names = []
        self.l_outvars = []
        self.d_outvars = {}
        self.l_connections = []

        self.gseim_prm = {}

#       assigning pad_in_names, pad_out_names: we need to
#       open the grc file for the subckt
#
#       category: '[GRC Hier Blocks]'
#       indicates subckt

        with open(os.path.expanduser(filename)) as f:
            data = yaml.load(f, Loader=yaml.FullLoader)

        self.flag_subckt = False
        self.flag_block = False

        if 'category' in data.keys():
            if 'Hier' in data['category']:
                self.flag_subckt = True
            else:
                self.flag_block = True
        else:
            self.flag_block = True

        if not data['id'].startswith('pad'):
            if 'parameters' in data.keys():
                for d in data['parameters']:
                    self.gseim_prm[d['id']] = d['default']

        if 'inputs' in data.keys():
            self.n_in = len(data['inputs'])
        else:
            self.n_in = 0
        if 'outputs' in data.keys():
            self.n_out = len(data['outputs'])
        else:
            self.n_out = 0

        if self.flag_subckt:
#           assign pad_in_names, pad_out_names
            grc_filename = data['grc_source']

            with open(os.path.expanduser(grc_filename)) as f:
                data = yaml.load(f, Loader=yaml.FullLoader)
            data = treat_virtual(data)
            for i in data['blocks']:
                if 'name' in i.keys():
                    if i['name'].startswith('pad_source'):
                        self.pad_in_names.append(i['name'])
                    if i['name'].startswith('pad_sink'):
                        self.pad_out_names.append(i['name'])

            self.l_connections = data['connections']
            self.d_outvars = data['outvars']

        if 'outvars' in data.keys():
            self.l_outvars = data['outvars']

class cct:
    def __init__(self, name, flag_top):
        self.name = name
        self.netlist_name = ''
        self.flag_top = False
        self.flag_subckt = False
        self.flag_block = True
        self.children = []

        if (flag_top):
            self.flag_top = True
            self.pos = []
            self.pos.append(0)
            self.pos.append(0)
        else:
            self.pos = []
        self.parent = None
        self.lib_map = -1
        self.l_in_names = []
        self.l_out_names = []
        self.l_connections = []
        self.d_outvars = {}
        self.index_sub = -1
        self.flag_prm_file = False
        self.prm_file_name = ''
        self.grc_file_name = ''

        self.node_map = []
        self.gseim_prm = {}

    def add_child(self, c):
        c.pos.append(self.pos[0] + 1)
        c.pos.append(len(self.children))
        self.children.append(c)
        c.parent = self
        self.flag_block = False

def parse_1(parent, child_name, dir_block, dir_sub, n_sub):
#
#   yyy: name of child (after removing $xxx)
#   dir_block: directory with all blocks (elements).
#   dir_sub: directory with all subckt's
#

    name1 = child_name.split('$')[0]
    filename = dir_block + name1 + '.block.yml'

    if (os.path.exists(filename)):
        cct0 = cct(child_name, False)
        parent.add_child(cct0)
        return
    else:
        filename = dir_sub + name1 + '.hblock.yml'
        if (os.path.exists(filename)):
            cct0 = cct(child_name, False)
            parent.add_child(cct0)
            cct0.flag_subckt = True
            cct0.index_sub = n_sub[0]
            n_sub[0] += 1
            if cct0.pos[0] > n_sub[1]:
                n_sub[1] = cct0.pos[0]
#           get elements called by this subckt:
            with open(os.path.expanduser(filename)) as f:
                data = yaml.load(f, Loader=yaml.FullLoader)

            grc_filename = data['grc_source']
            cct0.grc_file_name = grc_filename

            with open(os.path.expanduser(grc_filename)) as f:
                data = yaml.load(f, Loader=yaml.FullLoader)
            data = treat_virtual(data)
            l_connections = data['connections']
            l_unique = get_unique_names_1(l_connections)
            for i in l_unique:
                parse_1(cct0, i, dir_block, dir_sub, n_sub)
        else:
            print('parse_1:', child_name, 'is not block or subckt.', flush=True)
            print('   Halting...', flush=True)
            sys.exit()

def make_lib(input_entity, l_lib, dir_block, dir_sub):
#
#   We will assume that block names and subckt names are different.

    name1 = input_entity.name.split('$')[0]
    l_names = list(map(lambda x : x.name, l_lib))

    if input_entity.flag_subckt:
        if name1 in l_names:
            map0 = l_names.index(name1)
            input_entity.lib_map = map0
        else:
            filename = dir_sub + name1 + '.hblock.yml'
            lib0 = lib_entity(name1, filename)
            map0 = len(l_lib)
            l_lib.append(lib0)
            input_entity.lib_map = map0

        n_nodes = len(l_lib[map0].l_connections)
        for i in range(n_nodes):
            input_entity.node_map.append(None)
        for i in input_entity.children:
            make_lib(i, l_lib, dir_block, dir_sub)
    else:
        if name1 in l_names:
            map0 = l_names.index(name1)
        else:
            filename = dir_block + name1 + '.block.yml'
            lib0 = lib_entity(name1, filename)
            map0 = len(l_lib)
            l_lib.append(lib0)
        input_entity.lib_map = map0

def make_dict_1(input_entity, d1):
    if input_entity.flag_subckt:
        n_next = max(y for (x,y) in d1.items()) + 1
        d1[input_entity.name] = n_next
        for j in input_entity.children:
            make_dict_1(j, d1)
    else:
        return

def add_in_out(l_connections, l_c):
    for c in l_connections:
        l_c1 = [c[0], c[1], 'out']
        l_c2 = [c[2], c[3], 'in']
        l = [l_c1, l_c2]
        l_c.append(l)

def get_unique_wires(l_c, wires):
    wires.append(l_c[0])

    for i in range(1, len(l_c)):
        l1 = l_c[i]
        flag_found = False
        for i_w in range(0, len(wires)):
#           check if the first element is in the wire
            if l1[0] in wires[i_w]:
#               add the second element to the wire
                wires[i_w].append(l1[1])
                flag_found = True
                break
#           check if the second element is in the wire
            if l1[1] in wires[i_w]:
#               add the first element to the wire
                wires[i_w].append(l1[0])
                flag_found = True
                break
#           check if the first element has a dummy source
            if l1[0][0].startswith('dummy_source'):
                wire_sink = [l1[0][0].replace('source', 'sink'), '0', 'in']
                if wire_sink in wires[i_w]:
                    wires[i_w].append(l1[0])
                    wires[i_w].append(l1[1])
                    flag_found = True
                    break
#           check if the second element has a dummy sink
            if l1[1][0].startswith('dummy_sink'):
                wire_source = [l1[1][0].replace('sink', 'source'), '0', 'out']
                if wire_source in wires[i_w]:
                    wires[i_w].append(l1[0])
                    wires[i_w].append(l1[1])
                    flag_found = True
                    break

        if flag_found:
            continue
        else:
#           new wire:
            wires.append(l1)

def assign_node_names(input_entity, l_lib, cct_file,
    dir_block, dir_sub, d_names, node_counter):
    if input_entity.flag_top:
        l_connections = input_entity.l_connections
    else:
        l_connections = l_lib[input_entity.lib_map].l_connections

    l_c = []
    add_in_out(l_connections, l_c)

    wires = []
    get_unique_wires(l_c, wires)

    for i in input_entity.children:
        for i_in in range(l_lib[i.lib_map].n_in):
            node0 = [i.name, str(i_in), 'in']
            flag_found = False
            for i_wire, wire in enumerate(wires):
                if node0 in wire:
                    flag_found = True
                    node_name = 'n' \
                      + str(d_names[i.parent.name]) \
                      + '_' + str(i_wire)
                    i.l_in_names.append(node_name)
            if not flag_found:
                i.l_in_names.append('n' + str(node_counter))
                node_counter += 1

        for i_out in range(l_lib[i.lib_map].n_out):
            node0 = [i.name, str(i_out), 'out']
            flag_found = False
            for i_wire, wire in enumerate(wires):
                if node0 in wire:
                    flag_found = True
                    node_name = 'n' \
                      + str(d_names[i.parent.name]) \
                      + '_' + str(i_wire)
                    i.l_out_names.append(node_name)
            if not flag_found:
                i.l_out_names.append('n' + str(node_counter))
                node_counter += 1

    for i in input_entity.children:
        if i.flag_subckt:
            assign_node_names(i, l_lib, cct_file,
               dir_block, dir_sub, d_names, node_counter)

def assign_gseim_prm(input_entity, l_lib, cct_file, dir_block, dir_sub):
    if input_entity.flag_top:
#       main cct: take parms from grc file of the cct:
        filename = cct_file
        with open(os.path.expanduser(filename)) as f:
            data = yaml.load(f, Loader=yaml.FullLoader)
        input_entity.gseim_prm = data['gparms']
    else:
#       subckt: take parms from grc file of the parent:
        lib_name = l_lib[input_entity.lib_map].name
        filename = input_entity.parent.grc_file_name
        with open(os.path.expanduser(filename)) as f:
            data = yaml.load(f, Loader=yaml.FullLoader)

        for d in data['blocks']:
            if d['name'] == input_entity.name:
                if 'parameters' in d.keys():
                    for k, v in d['parameters'].items():
                        if k in l_lib[input_entity.lib_map].gseim_prm.keys():
                            input_entity.gseim_prm[k] = v

    if input_entity.flag_subckt or input_entity.flag_top:
        for i in input_entity.children:
            assign_gseim_prm(i, l_lib, cct_file, dir_block, dir_sub)

def get_pad_in_name(subckt, input_pad_name, l_lib):
#   subckt is of class cct (it is the subckt which has input_pad_name
#   as one of its input pads).
#   input_pad_name is the name of the input pad, node_name (output)
#   is the node name of the pad in the subckt.

    index1 = l_lib[subckt.lib_map].pad_out_names.index(input_pad_name)
    node_name = subckt.l_out_names[index1]
    return node_name

def get_pad_out_name(subckt, output_pad_name, l_lib):
    index1 = l_lib[subckt.lib_map].pad_in_names.index(output_pad_name)
    node_name = subckt.l_in_names[index1]
    return node_name

def replace_node_names_1(input_entity, node_name, name_to_replace):
    for i in input_entity.parent.children:
        i.l_in_names = [node_name if x==name_to_replace else x
           for x in i.l_in_names]
        i.l_out_names = [node_name if x==name_to_replace else x
           for x in i.l_out_names]

def replace_nodes_1(input_entity, l_lib):
    if input_entity.flag_subckt:
        for i in input_entity.children:
            replace_nodes_1(i, l_lib)
    else:
        if input_entity.name.startswith('pad_source'):
#           l_out_names[0] since a pad has only one node
            node_name = get_pad_out_name(input_entity.parent,
               input_entity.name, l_lib)
            name_to_replace = input_entity.l_out_names[0]
            replace_node_names_1(input_entity, node_name, name_to_replace)
        elif input_entity.name.startswith('pad_sink'):
            node_name = get_pad_in_name(input_entity.parent,
               input_entity.name, l_lib)
            name_to_replace = input_entity.l_in_names[0]
            replace_node_names_1(input_entity, node_name, name_to_replace)

def make_netlist_name(input_entity, l_lib):
    if input_entity.flag_subckt:
        for i in input_entity.children:
            make_netlist_name(i, l_lib)
    else:
        if not input_entity.name.startswith('pad'):
            if input_entity.parent.flag_subckt:
                name1 = 's' + str(input_entity.parent.index_sub) \
                         + '#' + input_entity.name
            else:
                name1 = input_entity.name
            input_entity.netlist_name = name1

def make_gseim_line(input_entity, l_lib, l_lines):
#   l_lines is the list of lists, each entry corresponds to a line in
#   the gseim circuit file, e.g.,
#   [
#    ['xelement', 'name=xxx', 'type=xxx', 'x1=xx', ..],
#    ['xelement', 'name=xxx', 'type=xxx', 'x1=xx', ..],
#   ]

    if input_entity.flag_subckt:
        for i in input_entity.children:
            make_gseim_line(i, l_lib, l_lines)
    else:
        l_line = ['xelement']

        if not input_entity.name.startswith('pad'):
            name1 = input_entity.netlist_name

            l_line.append('name=' + name1)
            l_line.append('type=' + input_entity.name.split('$')[0])

            for i_node, node_name in enumerate(input_entity.l_in_names):
                lib_node_name = l_lib[input_entity.lib_map].gseim_in_names[i_node]
                l_line.append(lib_node_name + '=' + node_name)
            for i_node, node_name in enumerate(input_entity.l_out_names):
                lib_node_name = l_lib[input_entity.lib_map].gseim_out_names[i_node]
                l_line.append(lib_node_name + '=' + node_name)
            for k, v in input_entity.gseim_prm.items():
#               We will reserve keyword 'name' for GUI use only.
#               Skip it here.
                v1 = v.strip()
                if k != 'name':
                    if v1 != l_lib[input_entity.lib_map].gseim_prm[k]:
                        l_line.append(k + '=' + v1)

            l_lines.append(l_line)

def make_prm_files(input_entity, l_lib, dir_block, dir_sub, l_files):
    if input_entity.flag_subckt:
        name1 = input_entity.name.split('$')[0]
        filename = dir_sub + name1 + '_parm.py'
        if os.path.exists(filename):
            input_entity.flag_prm_file = True
            input_entity.prm_file_name = filename
            l_files.append(filename)

        for i in input_entity.children:
            make_prm_files(i, l_lib, dir_block, dir_sub, l_files)
    else:
        return

def compute_prm_0(f, d_prm):
    import compute_prm as comp
    method_to_call = getattr(comp, f)
    method_to_call(d_prm)

def compute_prm_1(input_entity, d_prm):
    if input_entity.flag_top:
        name1 = input_entity.prm_file_name.split('.')[0]
        f = name1.split('/')[-1]
    else:
        name1 = input_entity.name.split('$')[0]
        f = name1.split('/')[-1] + '_parm'

#   need to remove path from name1, keeping only the last part:
    compute_prm_0(f, d_prm)

    for prm in d_prm.values():
        if not isinstance(prm, str):
            print('compute_prm_1: prm is expected to be type str', flush=True)
            print('  check this:' + prm, flush=True)
            sys.exit()

def compute_prm_2(input_entity):
    if input_entity.flag_subckt:
#       assign parameters coming from parent:
        for k, v in input_entity.gseim_prm.items():
            if v in input_entity.parent.gseim_prm.keys():
                input_entity.gseim_prm[k] = input_entity.parent.gseim_prm[v]

#       evaluate parameters for this entity:
        if input_entity.flag_prm_file:
            compute_prm_1(input_entity, input_entity.gseim_prm)

#       check if any parameter = 'computed'
        if 'computed' in input_entity.gseim_prm.values():
            print('compute_prm_2: some parameter has not been computed.',
              'Halting...', flush=True)
            sys.exit()

#       treat children:
        for i in input_entity.children:
            compute_prm_2(i)
    else:
#       assign parameters coming from parent:
        for k, v in input_entity.gseim_prm.items():
            if v in input_entity.parent.gseim_prm.keys():
                input_entity.gseim_prm[k] = input_entity.parent.gseim_prm[v]

def assign_node_map(input_entity, l_lib):
    if input_entity.flag_top:
        l_connections = input_entity.l_connections
    else:
        l_connections = l_lib[input_entity.lib_map].l_connections

    if input_entity.flag_top or input_entity.flag_subckt:
        for i_c, connection in enumerate(l_connections):
            name1 = connection[0]
            in_node = int(connection[1])

            for i in input_entity.children:
                if i.name == name1:
                    input_entity.node_map[i_c] = i.l_out_names[in_node]

        for i in input_entity.children:
            if i.flag_subckt:
                assign_node_map(i, l_lib)

def resolve_outvar(ov_main, ov_temp, input_entity, d_ov, l_lib,
    d_flag_resolved, d_resolved_outvars):

    ov_type = d_ov[ov_temp][0]

    if ov_type == 'connection':
        connection = d_ov[ov_temp][1]
        if input_entity.flag_top:
            index1 = input_entity.l_connections.index(connection)
        elif input_entity.flag_subckt:
            index1 = l_lib[input_entity.lib_map].l_connections.index(connection)
        else:
            print('resolve_outvar: connection not resolved. Halting...',
              flush=True)
            sys.exit()

        d_flag_resolved[ov_main] = True
        d_resolved_outvars[ov_main] = 'xvar_of_' + input_entity.node_map[index1]
        return
    elif ov_type == 'outvar':
        e_name, ov_name = d_ov[ov_temp][1]

        if input_entity.flag_top:
            for child in input_entity.children:
                if child.name == e_name:
                    if child.flag_block:
                        d_flag_resolved[ov_main] = True
                        d_resolved_outvars[ov_main] = ov_name + '_of_' + e_name
                        return
                    else:
                        d0 = l_lib[child.lib_map].d_outvars
                        resolve_outvar(ov_main, ov_name, child, d0, l_lib,
                           d_flag_resolved, d_resolved_outvars)
        elif input_entity.flag_subckt:
            for child in input_entity.children:
                if child.name == e_name:
                    if child.flag_block:
                        d_flag_resolved[ov_main] = True
                        d_resolved_outvars[ov_main] = ov_name + '_of_' + child.netlist_name
                        return
                    else:
                        d0 = l_lib[child.lib_map].d_outvars
                        resolve_outvar(ov_main, ov_name, child, d0, l_lib,
                           d_flag_resolved, d_resolved_outvars)
        else:
            print('resolve_outvar: outvar not resolved (2). Halting...',
              flush=True)
            sys.exit()

# main program:

if len(sys.argv) != 6:
    print('gseim1.py: need 5 arguments. Halting...')
    sys.exit()

gseim_file = sys.argv[1]
dir_block  = sys.argv[2]
dir_sub    = sys.argv[3]
dir_xbe    = sys.argv[4]
cct_file   = sys.argv[5]

dir_exec = dir_xbe.replace('xbe', 'exec')

print('main: cct_file:', cct_file, flush=True)
with open(os.path.expanduser(cct_file)) as f:
    data = yaml.load(f, Loader=yaml.FullLoader)
data = treat_virtual(data)

cct1 = cct('top_cct', True)
cct1.grc_file_name = cct_file

cct1.l_connections = data['connections']
cct1.d_outvars = data['outvars']
l_top_elements = get_unique_names_1(cct1.l_connections)

n_main_nodes = len(cct1.l_connections)
for i in range(n_main_nodes):
    cct1.node_map.append(None)

# n_sub[0]: number of subckts in the circuit
# n_sub[1]: depth of the circuit
n_sub = [0, 0]
for i in l_top_elements:
    parse_1(cct1, i, dir_block, dir_sub, n_sub)

# prepare list of library entities (which includes basic
# elements and subcircuits)

l_lib = []
for i in cct1.children:
    make_lib(i, l_lib, dir_block, dir_sub)

# prepare dict (to be used in circuit node names):
d_names = {}
d_names[cct1.name] = 0
for i in cct1.children:
    make_dict_1(i, d_names)

node_counter = 0
assign_node_names(cct1, l_lib, cct_file,
    dir_block, dir_sub, d_names, node_counter)

# extract input/output nodes names from GSEIM templates.
# Assume the template name to be xxx.xbe where xxx comes
# from l_lib.

for i in l_lib:
    if not i.flag_subckt:
        name1 = i.name
        if not name1.startswith('pad'):
            filename = dir_xbe + name1 + '.xbe'
            gu.extract_strings_2(filename, 'input_vars:', i.gseim_in_names)
            gu.extract_strings_2(filename, 'output_vars:', i.gseim_out_names)
            gu.extract_dict_1(filename, 'rparms:', i.gseim_prm)

# replace node names (subckt parsing):
for depth in range(n_sub[1]):
    for i in cct1.children:
        replace_nodes_1(i, l_lib)

# make a list of files to be called for prm computation:
l_prm_files = []

# main circuit file:
filename = cct_file.replace('.grc', '_parm.py')
if os.path.exists(filename):
    cct1.flag_prm_file = True
    cct1.prm_file_name = filename
    l_prm_files.append(filename)

for i in cct1.children:
    make_prm_files(i, l_lib, dir_block, dir_sub, l_prm_files)

# concatenate the prm computation files:
# for subckt, s_1, expect file s_1_parm.py

filename = dir_exec + 'compute_prm.py'
with open(os.path.expanduser(filename), 'w') as f_out:
    for infile in l_prm_files:
        with open(os.path.expanduser(infile), 'r') as f_in:
            f_out.write(f_in.read())

assign_gseim_prm(cct1, l_lib, cct_file, dir_block, dir_sub)

# evaluate main cct prm expression, if file exisits:

if cct1.flag_prm_file:
    compute_prm_1(cct1, cct1.gseim_prm)

if 'computed' in cct1.gseim_prm.values():
    print('main: some parameter has not been computed. Halting...', flush=True)
    sys.exit()

# evaluate prm expression for rest of the tree:

for i in cct1.children:
    compute_prm_2(i)

assign_node_map(cct1, l_lib)

# prepare names of leaf blocks as they would appear in the netlist
# (and also in outvars)

for i in cct1.children:
    make_netlist_name(i, l_lib)

n_outvars = len(cct1.d_outvars)
d_resolved_outvars = {}
d_flag_resolved = {}
for k in cct1.d_outvars.keys():
    d_flag_resolved[k] = False

for k in cct1.d_outvars.keys():
    resolve_outvar(k, k, cct1, cct1.d_outvars, l_lib,
       d_flag_resolved, d_resolved_outvars)

for k in cct1.d_outvars.keys():
    if not d_flag_resolved[k]:
        print('main: outvar', k, ' not resolved. Halting...', flush=True)
        sys.exit()

# prepare list of lines to be written to gseim file:
l_lines = []

for i in cct1.children:
    make_gseim_line(i, l_lib, l_lines)

d_slvlib = {}
filename = dir_exec + 'gseim_slvparms.in'
gu.get_parms_1(filename, d_slvlib)
l_solve_blocks = sorted(data['solve_blocks'], key=lambda x: int(x['index']))
l_output_blocks = data['output_blocks']

nmax = 80
indent1 = '   '
indent2 = '+     '
indent3 = '     '
indent4 = '+       '
f1 = open(os.path.expanduser(gseim_file), 'w')
f1.write('title: ' + cct_file + '\n')
f1.write('begin_circuit\n')

for line in l_lines:
    if not 'dummy_' in line[2]:
        gu.format_string_4a(f1, nmax, indent1, indent2, line)

for k,v in d_resolved_outvars.items():
    f1.write(indent1 + 'outvar: ' + k + '=' + v + '\n')

f1.write('end_circuit\n')

for slv in l_solve_blocks:
    f1.write('begin_solve\n')
    f1.write(indent1 + 'solve_type=' + slv['d_parms']['solve_type'] + '\n')
    f1.write(indent1 + 'initial_sol=' + slv['d_parms']['initial_sol'] + '\n')
    if slv['d_parms']['initial_sol'] == 'read_from_file':
        f1.write(indent1 + 'initial_sol_file='
          + slv['d_parms']['initial_sol_file'] + '\n')

    l_skip = [
      'solve_type',
      'initial_sol',
      'initial_sol_file',
      'output_solution_file',
      'block_index',
    ]

    for k, v in slv['d_parms'].items():
        if k not in l_skip:
            l2 = [d_slvlib[k]['default'], 'none']
            if v in l2:
                if d_slvlib[k]['flag_force_write']:
                    v1 = cct1.gseim_prm[v] if v in cct1.gseim_prm.keys() else v
                    f1.write(indent1 + 'method: ' + k + '=' + v1 + '\n')
            else:
                v1 = cct1.gseim_prm[v] if v in cct1.gseim_prm.keys() else v
                f1.write(indent1 + 'method: ' + k + '=' + v1 + '\n')

    for out_name in slv['output_blocks']:
        l_out_1 = list(filter(lambda x: x['name'] == out_name, l_output_blocks))
        if len(l_out_1) != 1:
            print('did not find output block', out_name,
              'in circuit file. Halting...', flush=True)
            sys.exit()
        out = l_out_1[0]
        f1.write(indent1 + 'begin_output\n')
        if out['d_parms']['filename'] == 'none':
            print('filename not specified in output block',
              out.name, 'Halting...', flush=True)
            sys.exit()

        l_line = []
        l_line.append('filename=' + out['d_parms']['filename'])
        if out['d_parms']['append'] == 'yes':
            l_line.append('append=yes')

        l_line.append('limit_lines=' + out['d_parms']['limit_lines'])

        gu.format_string_4a(f1, nmax, indent3, indent4, l_line)

        l_control = ['control:']
        if out['d_parms']['fixed_interval'] != 'none':
            l_control.append('fixed_interval=' + out['d_parms']['fixed_interval'])
        if out['d_parms']['t_start'] != 'none':
            l_control.append('out_tstart=' + out['d_parms']['t_start'])
        if out['d_parms']['t_end'] != 'none':
            l_control.append('out_tend=' + out['d_parms']['t_end'])

        if len(l_control) > 1:
            gu.format_string_4a(f1, nmax, indent3, indent4, l_control)

        if not out['outvars']:
            print('no outvars specified in output block',
              out.name, 'Halting...', flush=True)
            sys.exit()
        l_line = []
        l_line.append('variables:')
        for ov in out['outvars']:
            l_line.append(ov)
        gu.format_string_4a(f1, nmax, indent3, indent4, l_line)

        f1.write(indent1 + 'end_output\n')

    if slv['d_parms']['output_solution_file'] != 'none':
        f1.write(indent1 + 'begin_output\n')
        f1.write(indent3 + 'filename=' + slv['d_parms']['output_solution_file'] + '\n')
        f1.write(indent3 + 'variables: solution' + '\n')
        f1.write(indent1 + 'end_output\n')

    f1.write('end_solve\n')

f1.write('end_cf\n')
f1.close()

# check if there are filter_1 elements in the file

f1 = open(os.path.expanduser(gseim_file), 'r')
flag_filter_1 = False
while True:
    line = f1.readline()
    l = line.split()
    if l[0] == 'end_cf':
        break
    elif l[0] == 'xelement':
        if l[1].split('=')[-1].startswith('filter_1'):
            flag_filter_1 = True
            break
f1.close()
if flag_filter_1:
    cmd = 'python3 ' + dir_exec + 'parse_filters.py ' + gseim_file + ' ' + dir_exec
    print('filter_1 found: cmd:', cmd, flush=True)
    os.system(cmd)
    cmd = 'mv cct_temp_1.in ' + gseim_file
    os.system(cmd)

print('Program completed.', flush=True) 
