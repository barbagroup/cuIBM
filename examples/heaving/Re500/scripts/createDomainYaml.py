"""
Creates a YAML file with info about the structured Cartesian mesh that will be
parsed by cuIBM.
"""

from snake.cartesianMesh import CartesianStructuredMesh


# info about the 2D structured Cartesian grid
info = [{'direction': 'x', 'start': -8.0,
         'subDomains': [{'end': -0.52,
                         'width': 0.005,
                         'stretchRatio': 1.02,
                         'reverse': True},
                        {'end': 0.52,
                         'width': 0.005,
                         'stretchRatio': 1.0},
                        {'end': 0.78,
                         'width': 0.005,
                         'stretchRatio': 1.02},
                        {'end': 10.0,
                         'width': 0.01,
                         'stretchRatio': 1.0}]},
        {'direction': 'y', 'start': -12.0,
         'subDomains': [{'end': -0.52,
                         'width': 0.005,
                         'stretchRatio': 1.015,
                         'reverse': True},
                        {'end': 0.52,
                         'width': 0.005,
                         'stretchRatio': 1.0},
                        {'end': 12.0,
                         'width': 0.005,
                         'stretchRatio': 1.015}]}]

mesh = CartesianStructuredMesh()
mesh.create(info, mode='cuibm')
mesh.print_parameters()
mesh.write_yaml_file('domain.yaml')
