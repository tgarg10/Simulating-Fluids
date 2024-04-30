from vpython import *
import numpy as np

import globals

class particles_setup:
    def __init__(self):
        self.particles_list = np.array([])
        self.densitites = np.array([])

    # Setting the valls in a grid formation
    def start_arrangement(self):

        # Place particles in a grid formation
        particles_per_row = int(np.sqrt(globals.total_particles))
        particles_per_col = int((globals.total_particles - 1) / particles_per_row) + 1

        for i in range(globals.total_particles):
            particle = sphere(color = color.cyan, radius = globals.particle_radius, make_trail=False, retain=200)
            particle.mass = globals.particle_mass
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

    # Find the influence of a particle at a distance from it.
    # Influence reduces with distance from the center of the particle
    def smoothing_kernel(self, distance):
        volume = np.pi * globals.influence_radius ** 8 / 4
        value = np.max(0, globals.influence_radius ** 2 - distance ** 2)
        return value ** 3 / volume

    def smoothing_kernel_derivtive(self, distance):
        if (distance > globals.influence_radius):
            return 0.0
        f = globals.influence_radius ** 2 - distance ** 2
        scale =- -24 / (np.pi * globals.influence_radius ** 8)
        return scale * distance * f ** 2

    # Calculate the density of at a specific point
    # to move them from areas of high density to low density
    def calculate_density(self, sample_point):
        density = 0

        # Loop over all particle positions to get the density at a point
        for particle in self.particles_list:
            distance = mag(particle.pos - sample_point)
            influence = self.smoothing_kernel(distance)
            density += globals.particle_mass * influence

        return density

    # Pre-calculting the values of densities
    def update_densities(self):
        for i in range(globals.total_particles):
            self.densitites[i] = self.calculate_density(self.particles_list[i])

        return self.densitites

    # Using gradient descent to get the particles to a zone of low density
    def calculate_pressure_force(self, sample_point):
        pressure_force = vector(0, 0, 0)
        self.update_densities()

        # Loop over all particle positions to get the density at a point
        for i in range(globals.total_particles):
            particle = self.particles_list[i]
            distance = mag(particle.pos - sample_point)
            direction = (particle.pos - sample_point) / distance
            slope = self.smoothing_kernel_derivative(distance)
            density = self.densities[i]
            pressure_force += -self.convert_density_to_pressure(density) direction * slope * globals.particle_mass / density

        return pressure_force
    
    # Pressure from a particle to another
    def convert_density_to_pressure(density):
        density_error = density - globals.target_density
        pressure = density_error * globals.pressure_multiplier
        return pressure