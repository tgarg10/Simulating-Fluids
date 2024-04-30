# Global Variables ---

# Setting up the wall around the Fluid
initial_wall_length = 10.0
initial_wall_height = 10.0
initial_wall_width = 8.0
wall_thickness = 0.1

# particle parameters
gravity = 9.8
update_rate = 200
delta_time = 0.01

collision_damping = 0.9
particle_radius = 0.2
total_particles = 200

# Spacing between two particles for the starting grid
spacing = 0.4

# Keeping the particles away from the wall
wall_length = initial_wall_length - wall_thickness * 0.5 - particle_radius
wall_height = initial_wall_height - wall_thickness * 0.5 - particle_radius
wall_width = initial_wall_width - wall_thickness * 0.5 - particle_radius