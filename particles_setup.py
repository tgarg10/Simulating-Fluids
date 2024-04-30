from vpython import *
import numpy as np

import globals

class particles_setup:
    def __init__(self):
        self.particles_list = np.array([])

    # Setting the valls in a grid formation
    def start_arrangement(self):

        # Place particles in a grid formation
        particles_per_row = int(np.sqrt(globals.total_particles))
        particles_per_col = int((globals.total_particles - 1) / particles_per_row) + 1

        for i in range(globals.total_particles):
            particle = sphere(color = color.cyan, radius = globals.particle_radius, make_trail=False, retain=200)
            particle.mass = 1.0
            particle.p = vector(0, 0, 0)
            particle.pos.x = (i % particles_per_row - particles_per_row / 2 + 0.5) * globals.spacing
            particle.pos.y = (i / particles_per_row - particles_per_col / 2 + 0.5) * globals.spacing
            self.particles_list = np.append(self.particles_list, particle)
        
        return self.particles_list

    # Moving the particles
    def move_particles(self):
        # Iterating through the list of particles
        for i in range(len(self.particles_list)):
            self.add_gravity(self.particles_list[i])
            self.resolve_collisions(self.particles_list[i], globals.collision_damping)


    # Adding gravity to the particle's movement
    def add_gravity(self, particle):
        particle.p.y += -1 * globals.gravity * globals.delta_time
        particle.pos.y += particle.p.y * globals.delta_time


    # Reflects the particle off if the particle collides the box
    def resolve_collisions(self, particle, collision_damping):
        if not (globals.wall_length > particle.pos.x > -globals.wall_length):
            particle.pos.x = globals.wall_length * np.sign(particle.pos.x)
            particle.p.x *= -1 * collision_damping
        if not (globals.wall_height > particle.pos.y > -globals.wall_height):
            particle.pos.y = globals.wall_height * np.sign(particle.pos.y)
            particle.p.y *= -1 * collision_damping
        if not (globals.wall_width > particle.pos.z > -globals.wall_width):
            particle.pos.z = globals.wall_width * np.sign(particle.pos.z)
            particle.p.z *= -1 * collision_damping
        