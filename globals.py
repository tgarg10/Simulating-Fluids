# Global Variables ---

# Setting up the wall around the Fluid
initial_wall_length = 4.0
initial_wall_height = 4.0
initial_wall_width = 4.0
wall_thickness = 0.1

# particle parameters
gravity = 9.8
update_rate = 200
delta_time = 0.1
particle_mass = 1.0

collision_damping = 0.9
particle_radius = 0.5
influence_radius = particle_radius * 4 # Radius of particle's Influence
total_particles = 100

# Spacing between two particles for the starting grid
spacing = 0.4

target_density = 2.75
pressure_multiplier = 10


# Keeping the particles away from the wall
wall_length = initial_wall_length - wall_thickness * 0.5 - particle_radius
wall_height = initial_wall_height - wall_thickness * 0.5 - particle_radius
wall_width = initial_wall_width - wall_thickness * 0.5 - particle_radius
