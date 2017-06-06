"""
From a boundary (2D closed line), uniformly discretizes the interior with a
given grid-spacing.
"""

import os

from snake.geometry import Geometry


filepath = os.path.join(os.environ['SNAKE'],
                        'resources',
                        'geometries',
                        'flyingSnake2d',
                        'flyingSnake2dAoA35.dat')
body = Geometry(file_path=filepath, skiprows=1)
body.keep_inside(ds=0.004)
body.write(file_path='new_body.txt')
