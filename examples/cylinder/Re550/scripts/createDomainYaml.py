"""
Creates a YAML file with info about the structured Cartesian mesh that will be
parsed by cuIBM.
"""

from snake.cartesianMesh import CartesianStructuredMesh


# info about the 2D structured Cartesian grid
width = 0.01  # minimum grid spacing in the x- and y- directions
info = [{'direction': 'x', 'start': -15.0,
         'subDomains': [{'end': -0.54,
                         'width': width,
                         'stretchRatio': 1.02,
                         'reverse': True},
                        {'end': 0.54,
                         'width': width,
                         'stretchRatio': 1.0},
                        {'end': 15.0,
                         'width': width,
                         'stretchRatio': 1.02}]},
        {'direction': 'y', 'start': -15.0,
         'subDomains': [{'end': -0.54,
                         'width': width,
                         'stretchRatio': 1.02,
                         'reverse': True},
                        {'end': 0.54,
                         'width': width,
                         'stretchRatio': 1.0},
                        {'end': 15.0,
                         'width': width,
                         'stretchRatio': 1.02}]}]

mesh = CartesianStructuredMesh()
mesh.create(info, mode='cuibm')
mesh.print_parameters()
mesh.write_yaml_file('domain.yaml')
