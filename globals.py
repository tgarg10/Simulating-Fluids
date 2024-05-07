# Global Variables ---

# Setting up the wall around the Fluid
initial_wall_length = 3.0
initial_wall_height = 2.0
initial_wall_width = 0.4
wall_thickness = 0.1

# particle parameters
gravity = 4
update_rate = 1000
delta_time = 0.01
particle_mass = 1.0

collision_damping = 0.0
particle_radius = 0.3
influence_radius = 1.2 * particle_radius # Radius of particle's Influence
total_particles = 40

# Spacing between two particles for the starting grid
spacing = 0.2

target_density = 1.75
viscosity_strength = 0.05
pressure_multiplier = 40


# Keeping the particles away from the wall
wall_length = initial_wall_length - wall_thickness * 0.5 - particle_radius
wall_height = initial_wall_height - wall_thickness * 0.5 - particle_radius
wall_width = initial_wall_width - wall_thickness * 0.5 - particle_radius


# Keeping the particles away from the wall
wall_length = initial_wall_length - wall_thickness * 0.5 - particle_radius
wall_height = initial_wall_height - wall_thickness * 0.5 - particle_radius
wall_width = initial_wall_width - wall_thickness * 0.5 - particle_radius