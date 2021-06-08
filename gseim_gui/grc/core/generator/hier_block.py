import collections
import os

import six
import codecs

from .top_block import TopBlockGenerator

from .. import Constants
from ..io import yaml

class HierBlockGenerator(TopBlockGenerator):
    """Extends the top block generator to also generate a block YML file"""

    def __init__(self, flow_graph, file_path):
        """
        Initialize the hier block generator object.

        Args:
            flow_graph: the flow graph object
            file_path: where to write the py file (the yml goes into HIER_BLOCK_LIB_DIR)
        """
        TopBlockGenerator.__init__(self, flow_graph, file_path)
        platform = flow_graph.parent

        hier_block_lib_dir = platform.config.hier_block_lib_dir
        if not os.path.exists(hier_block_lib_dir):
            os.mkdir(hier_block_lib_dir)

        self._mode = Constants.HIER_BLOCK_FILE_MODE
        print('HierBlockGenerator: hier_block_lib_dir:', hier_block_lib_dir)
        self.file_path = os.path.join(hier_block_lib_dir, self._flow_graph.get_option('id'))
        self.file_path_yml = self.file_path + '.hblock.yml'
        print('HierBlockGenerator: self.file_path_yml:', self.file_path_yml)

    def write(self, flow_graph):
        """generate output and write it to files"""
        print('HierBlockGenerator: calling TopBlockGenerator.write')
        TopBlockGenerator.write(self)
        print('HierBlockGenerator: TopBlockGenerator.write over')

        data = yaml.dump(self._build_block_n_from_flow_graph_io(flow_graph))

        replace = [
            ('parameters:', '\nparameters:'),
            ('inputs:', '\ninputs:'),
            ('outputs:', '\noutputs:'),
            ('asserts:', '\nasserts:'),
            ('templates:', '\ntemplates:'),
            ('documentation:', '\ndocumentation:'),
        ]
        for r in replace:
            data = data.replace(*r)

        print('HierBlockGenerator: self.file_path_yml:', self.file_path_yml)
        with codecs.open(os.path.expanduser(self.file_path_yml), 'w', encoding='utf-8') as fp:
            fp.write(data)

        # Windows only supports S_IREAD and S_IWRITE, other flags are ignored
        os.chmod(self.file_path_yml, self._mode)

    def _build_block_n_from_flow_graph_io(self, flow_graph):
        """
        Generate a block YML nested data from the flow graph IO

        Returns:
            a yml node tree
        """
        # Extract info from the flow graph
        block_id = self._flow_graph.get_option('id')
        parameters = self._flow_graph.get_parameters()

        def var_or_value(name):
            if name in (p.name for p in parameters):
                return "${" + name + " }"
            return name

        # Build the nested data
        data = collections.OrderedDict()
        data['id'] = block_id
        data['label'] = (
            self._flow_graph.get_option('title') or
            self._flow_graph.get_option('id').replace('_', ' ').title()
        )
        data['category'] = self._flow_graph.get_option('category')

#       Parameters
        data['parameters'] = []

        if 'name' not in [dx['id'] for dx in data['parameters']]:
            p = collections.OrderedDict()
            p['id'] = 'name'
            p['label'] = 'name'
            p['dtype'] = 'string'
            p['default'] = 'none'
            data['parameters'].append(p)

        for k, v in flow_graph.gparms.items():
            p = collections.OrderedDict()
            p['id'] = k
            p['label'] = k
            p['default'] = v
            data['parameters'].append(p)

#       outvars
        data['outvars'] = []

        for k in flow_graph.outvars.keys():
            data['outvars'].append(k)

        # Ports
        for direction in ('inputs', 'outputs'):
            data[direction] = []
            for port in get_hier_block_io(self._flow_graph, direction):
                p = collections.OrderedDict()
                p['label'] = port.parent.params['label'].value
                if port.domain != Constants.DEFAULT_DOMAIN:
                    p['domain'] = port.domain
                p['dtype'] = port.dtype
                if port.domain != Constants.GR_MESSAGE_DOMAIN:
                    p['vlen'] = var_or_value(port.vlen)
                if port.optional:
                    p['optional'] = True
                data[direction].append(p)

        data['grc_source'] = str(self._flow_graph.grc_file_path)

        return data

def get_hier_block_io(flow_graph, direction, domain=None):
    """
    Get a list of io ports for this flow graph.

    Returns a list of blocks
    """
    pads = flow_graph.get_pad_sources() if direction == 'inputs' else flow_graph.get_pad_sinks()

    for pad in pads:
        for port in (pad.sources if direction == 'inputs' else pad.sinks):
            if domain and port.domain != domain:
                continue
            yield port
