from vpython import *

from particles_setup import particles_setup
from scene_setup import scene_setup
import globals

# Setting up things before the simulation starts
scene_setup()
particles_manager = particles_setup()
# Creating the list of particles
particles_list = particles_manager.start_arrangement()

# Update function
def update():
    particles_manager.move_particles()

while True:
    update()
    rate(globals.update_rate)    