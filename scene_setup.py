from vpython import *

import globals

class scene_setup:
    def __init__(self):
        # Instructions to move the scene.
        scene.caption = """Right button drag or Ctrl-drag to rotate "camera" to view scene.
        To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
            On a two-button mouse, middle is left + right.
        Shift-drag to pan left/right and up/down.
        Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""

        # Setting up the scene
        wallR = box(pos=vector(globals.initial_wall_length, 0, 0), size=vector(globals.wall_thickness, 2 * globals.initial_wall_height, 2 * globals.initial_wall_width),  color = color.white)
        wallL = box(pos=vector(-globals.initial_wall_length, 0, 0), size=vector(globals.wall_thickness, 2 * globals.initial_wall_height, 2 * globals.initial_wall_width),  color = color.white)
        wallB = box(pos=vector(0, -globals.initial_wall_height, 0), size=vector(2 * globals.initial_wall_length, globals.wall_thickness, 2 * globals.initial_wall_width),  color = color.white)
        wallT = box(pos=vector(0,  globals.initial_wall_height, 0), size=vector(2 * globals.initial_wall_length, globals.wall_thickness, 2 * globals.initial_wall_width),  color = color.white)
        wallBK = box(pos=vector(0, 0, -globals.initial_wall_width), size=vector(2 * globals.initial_wall_length, 2 * globals.initial_wall_height, globals.wall_thickness), color = color.white)
