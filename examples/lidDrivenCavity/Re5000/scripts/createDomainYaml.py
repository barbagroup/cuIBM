"""
Creates a YAML file with info about the structured Cartesian mesh that will be
parsed by cuIBM.
"""

from snake.cartesianMesh import CartesianStructuredMesh


# info about the 2D structured Cartesian grid
info = [{'direction': 'x', 'start': 0.0,
         'subDomains': [{'end': 1.0,
                         'width': 1.0 / 192,
                         'stretchRatio': 1.0}]},
        {'direction': 'y', 'start': 0.0,
         'subDomains': [{'end': 1.0,
                         'width': 1.0 / 192,
                         'stretchRatio': 1.0}]}]

mesh = CartesianStructuredMesh()
mesh.create(info, mode='cuibm')
mesh.print_parameters()
mesh.write_yaml_file('domain.yaml')
