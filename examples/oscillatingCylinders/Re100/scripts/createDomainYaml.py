"""
Creates a YAML file with info about the structured Cartesian mesh that will be
parsed by cuIBM.
"""

from snake.cartesianMesh import CartesianStructuredMesh


# info about the 2D structured Cartesian grid
info = [{'direction': 'x', 'start': -15.0,
         'subDomains': [{'end': -2.0,
                         'width': 0.02,
                         'stretchRatio': 1.02,
                         'reverse': True},
                        {'end': 2.0,
                         'width': 0.02,
                         'stretchRatio': 1.0},
                        {'end': 15.0,
                         'width': 0.02,
                         'stretchRatio': 1.02}]},
        {'direction': 'y', 'start': -15.0,
         'subDomains': [{'end': -1.0,
                         'width': 0.02,
                         'stretchRatio': 1.02,
                         'reverse': True},
                        {'end': 1.0,
                         'width': 0.02,
                         'stretchRatio': 1.0},
                        {'end': 15.0,
                         'width': 0.02,
                         'stretchRatio': 1.02}]}]

mesh = CartesianStructuredMesh()
mesh.create(info, mode='cuibm')
mesh.print_parameters()
mesh.write_yaml_file('domain.yaml')
