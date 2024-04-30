from vpython import *
import numpy as np

import globals

class particles_setup:
    def __init__(self):
        self.particles_list = np.array([])
        self.densities = np.zeros(globals.total_particles)
        self.predicted_positions = np.array([])
        self.spatial_lookup = dict()

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
            self.predicted_positions = np.append(self.predicted_positions, particle.pos)

        self.update_spatial_lookup()

        return self.particles_list

    # Moving the particles
    def move_particles(self):
        # Iterating through the list of particles

        # Apply gravity
        for i in range(globals.total_particles):
            self.particles_list[i].p += vector(0, -1 * globals.gravity * globals.delta_time, 0)
            # Predicting the particle's next position to better calculate its force
            self.predicted_positions[i] = self.particles_list[i].pos + self.particles_list[i].p * globals.delta_time

        self.update_spatial_lookup()

        # Calculate densities
        for i in range(globals.total_particles):
            self.densities[i] = self.calculate_density(self.predicted_positions[i])
        
        # Calculate and apply pressure forces
        for i in range(globals.total_particles):
            pressure_force = self.calculate_pressure_force(self.particles_list[i].pos)
            pressure_accelertion = pressure_force / self.densities[i]
            self.particles_list[i].p += pressure_accelertion * globals.delta_time

        # Update positions and resolve collisions
        for i in range(globals.total_particles): 
            self.particles_list[i].pos += self.particles_list[i].p * globals.delta_time           
            self.resolve_collisions(self.particles_list[i], globals.collision_damping)


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
        if (distance > globals.influence_radius):
            return 0.0
        volume = np.pi * globals.influence_radius ** 4 / 6
        return (globals.influence_radius - distance) ** 2 / volume


    def smoothing_kernel_derivtive(self, distance):
        if (distance > globals.influence_radius):
            return 0.0
        scale = 12 / (np.pi * globals.influence_radius ** 4)
        return (distance - globals.influence_radius) * scale


    # Calculate the density of at a specific point
    # to move them from areas of high density to low density
    def calculate_density(self, sample_point):
        density = 0

        # Loop over all particle positions to get the density at a point
        for particle in self.get_lookup_particles(sample_point):
            distance = mag(particle - sample_point)
            influence = self.smoothing_kernel(distance)
            density += globals.particle_mass * influence

        return density
    

    # Using gradient descent to get the particles to a zone of low density
    def calculate_pressure_force(self, sample_point):
        pressure_force = vector(0, 0, 0)
        # Loop over all particle positions to get the density at a point
        for i in range(globals.total_particles):
            particle = self.particles_list[i]
            distance = mag(particle.pos - sample_point)
            if distance != 0:
                direction = (particle.pos - sample_point) / distance
            else:
                direction = vector(random(), random(), random())
            slope = self.smoothing_kernel_derivtive(distance)
            density = self.densities[i]
            # Every action has an equal and opposite reaction
            shared_pressure = self.calculate_shared_pressure(density, self.densities[i])
            pressure_force += direction * shared_pressure * slope * globals.particle_mass / density

        return pressure_force
    
    # Assigning each particle to a cell in the spatial lookup
    def update_spatial_lookup(self):
        self.spatial_lookup = dict()
        for particle_pos in self.predicted_positions:
            particle_cell = self.position_to_cell_coord(particle_pos)
            if particle_cell in self.spatial_lookup:
                self.spatial_lookup[particle_cell] = np.append(self.spatial_lookup[particle_cell], particle_pos)
            else:
                self.spatial_lookup[particle_cell] = np.array([particle_pos])

        return self.spatial_lookup

    # Returns the 3x3 grid surrounding the particle
    def get_lookup_particles(self, particle_pos):
        cell_coord = self.position_to_cell_coord(particle_pos)
        influenced_particles = self.spatial_lookup[cell_coord]
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, 0, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, 1, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, 1, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, 0, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, 0, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, 1, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, 1, 1)))])

        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, 0, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, -1, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, -1, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, 0, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, 0, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, -1, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, -1, -1)))])

        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, 1, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (0, -1, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, 0, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, 0, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, -1, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, 1, 0)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, 1, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, -1, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, 1, 1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (1, -1, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, 1, -1)))])
        influenced_particles = np.concatenate([influenced_particles, self.check_and_add(self.vector_add(cell_coord, (-1, -1, 1)))])

        return influenced_particles

    # Check if the coordinate exists in the spatial lookup and return those particles
    def check_and_add(self, check_coord):
        if check_coord in self.spatial_lookup:
            return self.spatial_lookup[check_coord]
        else:
            return np.array([])

    # Pressure from a particle to another
    def convert_density_to_pressure(self, density):
        density_error = density - globals.target_density
        pressure = density_error * globals.pressure_multiplier
        return pressure

    # The average of the pressures of the two particles
    def calculate_shared_pressure(self, densityA, densityB):
        pressureA = self.convert_density_to_pressure(densityA)
        pressureB = self.convert_density_to_pressure(densityB)
        return (pressureA + pressureB) / 2
    
    # Returns the coordinate of the cell on the grid
    def position_to_cell_coord(self, particle_pos):
        offset_vector = vector(particle_pos.x % globals.influence_radius, 
                               particle_pos.y % globals.influence_radius, 
                               particle_pos.z % globals.influence_radius)
        cell_coord = (particle_pos - offset_vector) / globals.influence_radius
        return (cell_coord.x, cell_coord.y, cell_coord.z)
    
    def vector_add(self, p1, p2):
        return (p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2])