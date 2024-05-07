from vpython import *

# Global Variables ---

# Setting up the wall around the Fluid
initial_wall_length = 4.0
initial_wall_height = 3.0
initial_wall_width = 0.3
wall_thickness = 0.1

# particle parameters
gravity = 4
update_rate = 100
delta_time = 0.01
particle_mass = 1.0

collision_damping = 0.0
particle_radius = 0.2
influence_radius = 1.2 * particle_radius # Radius of particle's Influence
total_particles = 200

# Spacing between two particles for the starting grid
spacing = 0.2

target_density = 2.75
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

# Instructions to move the scene.
scene.caption = """Right button drag or Ctrl-drag to rotate "camera" to view scene.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
    On a two-button mouse, middle is left + right.
Shift-drag to pan left/right and up/down.
Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""

# Setting up the scene
wallR = box(pos=vector(initial_wall_length, 0, 0), size=vector(wall_thickness, 2 * initial_wall_height, 2 * initial_wall_width),  color = color.white)
wallL = box(pos=vector(-initial_wall_length, 0, 0), size=vector(wall_thickness, 2 * initial_wall_height, 2 * initial_wall_width),  color = color.white)
wallB = box(pos=vector(0, -initial_wall_height, 0), size=vector(2 * initial_wall_length, wall_thickness, 2 * initial_wall_width),  color = color.white)
wallT = box(pos=vector(0,  initial_wall_height, 0), size=vector(2 * initial_wall_length, wall_thickness, 2 * initial_wall_width),  color = color.white)
wallBK = box(pos=vector(0, 0, -initial_wall_width), size=vector(2 * initial_wall_length, 2 * initial_wall_height, wall_thickness), color = color.white)


particles_list = []
densities = [0.0 for i in range(total_particles)]
predicted_positions = []
spatial_lookup = dict()

# Setting the valls in a grid formation
def start_arrangement():
    global spatial_lookup
    # Place particles in a grid formation
    cube_width = floor(pow(total_particles, 1/3)) + 1

    for i in range(total_particles):
        particle = sphere(color = color.cyan, radius = particle_radius, make_trail=False, retain=200)
        particle.mass = particle_mass
        particle.p = vector(0, 0, 0)
        particle.pos.x = (i % cube_width) * spacing
        particle.pos.y = int(i / (cube_width ** 2)) * spacing
        particle.pos.z = (int(i / cube_width) % cube_width) * spacing
        particles_list.append(particle)
        predicted_positions.append(particle.pos)
    
    spatial_lookup = update_spatial_lookup()

    return particles_list

# Moving the particles
def move_particles():
    global spatial_lookup
    # Iterating through the list of particles

    # Apply gravity
    for i in range(total_particles):
        particles_list[i].p += vector(0, -1 * gravity * delta_time, 0)
        # Predicting the particle's next position to better calculate its force
        predicted_positions[i] = particles_list[i].pos + particles_list[i].p * delta_time

    spatial_lookup = update_spatial_lookup()

    # Calculate densities
    for i in range(total_particles):
        densities[i] = calculate_density(i)

    # Calculate and apply pressure forces
    for i in range(total_particles):
        pressure_force = calculate_pressure_force(i)
        pressure_force += calculate_viscosity_force(i)
        pressure_accelertion = pressure_force / densities[i]
        particles_list[i].p += pressure_accelertion * delta_time

    # Update positions and resolve collisions
    for i in range(total_particles): 
        particles_list[i].pos += particles_list[i].p * delta_time           
        resolve_collisions(particles_list[i], collision_damping)

# Reflects the particle off if the particle collides the box
def resolve_collisions(particle, collision_damping):
    if not (wall_length > particle.pos.x > -wall_length):
        particle.pos.x = wall_length * sign(particle.pos.x)
        particle.p.x *= -1 * collision_damping
    if not (wall_height > particle.pos.y > -wall_height):
        particle.pos.y = wall_height * sign(particle.pos.y)
        particle.p.y *= -1 * collision_damping
    if not (wall_width > particle.pos.z > -wall_width):
        particle.pos.z = wall_width * sign(particle.pos.z)
        particle.p.z *= -1 * collision_damping


# Find the influence of a particle at a distance from it.
# Influence reduces with distance from the center of the particle
def smoothing_kernel(distance):
    if (distance >= influence_radius):
        return 0.0
    volume = pi * influence_radius ** 4 / 6
    return (influence_radius - distance) ** 2 / volume

def smoothing_kernel_derivtive1(distance):
    if (distance >= influence_radius):
        return 0.0
    scale = 12 / (pi * influence_radius ** 4)
    return (distance - influence_radius) * scale

def smoothing_kernel_derivtive(distance):
    if (distance >= influence_radius):
        return 0.0
    f = influence_radius ** 2 - distance ** 2
    scale = -24 / (pi * influence_radius ** 8)
    return distance * scale * f ** 2

def viscosity_smoothing_kernel(distance):
    if (distance < influence_radius):
        volume = (64 * pi * influence_radius ** 9) / 315
        v = influence_radius ** 2 - distance ** 2
        return v ** 3 / volume
    return 0

# Calculate the density of at a specific point
# to move them from areas of high density to low density
def calculate_density(sample_index):
    density = 0

    # Loop over all particle positions to get the density at a point
    for other_index in get_lookup_particles(sample_index):
        distance = mag(particles_list[other_index].pos - particles_list[sample_index].pos)
        influence = smoothing_kernel(distance)
        density += particle_mass * influence

    return density

def calculate_viscosity_force(sample_index):
    viscosity_force = vector(0, 0, 0)

    # Loop over all particle positions to get the density at a point
    for other_index in get_lookup_particles(sample_index):
        distance = mag(particles_list[other_index].pos - particles_list[sample_index].pos)
        influence = viscosity_smoothing_kernel(distance)
        viscosity_force += (particles_list[other_index].p - particles_list[sample_index].p) * influence

    return viscosity_force * viscosity_strength

# Using gradient descent to get the particles to a zone of low density
def calculate_pressure_force(sample_index):
    pressure_force = vector(0, 0, 0)
    # Loop over all particle positions to get the density at a point
    for i in range(total_particles):
        if i == sample_index:
            continue
        particle = predicted_positions[i]
        distance = mag(particle - particles_list[sample_index].pos)
        if distance > 0.2:
            direction = (particle - particles_list[sample_index].pos) / distance
        else:
            direction = vector.random()
        slope = smoothing_kernel_derivtive(distance)
        density = densities[i]
        # Every action has an equal and opposite reaction
        shared_pressure = calculate_shared_pressure(density, densities[i])
        pressure_force += direction * shared_pressure * slope * particle_mass / density

    return pressure_force

# Assigning each particle to a cell in the spatial lookup
def update_spatial_lookup():
    spatial_lookup = dict()
    for i in range(total_particles):
        particle_cell = position_to_cell_coord(particles_list[i].pos)
        if particle_cell in spatial_lookup:
            spatial_lookup[particle_cell].append(i)
        else:
            spatial_lookup[particle_cell] = [i]
            
    return spatial_lookup

# Returns the 3x3 grid surrounding the particle
def get_lookup_particles(particle_index):
    cell_coord = position_to_cell_coord(particles_list[particle_index].pos)
    influenced_particles = spatial_lookup[cell_coord]

    offset_list = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (0, 0, -1), 
                    (0, -1, 0), (0, -1, -1), (-1, 0, 0), (-1, 0, -1), (-1, -1, 0), (-1, -1, -1), (0, 1, -1), (0, -1, 1), (1, 0, -1), 
                    (-1, 0, 1), (1, -1, 0), (-1, 1, 0), (1, 1, -1), (1, -1, 1), (-1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]

    for offset in offset_list:
        influenced_particles = influenced_particles + check_and_add(vector_add(cell_coord, offset))
    
    return influenced_particles

# Check if the coordinate exists in the spatial lookup and return those particles
def check_and_add(check_coord):
    if check_coord in spatial_lookup:
        return spatial_lookup[check_coord]
    else:
        return []

# Pressure from a particle to another
def convert_density_to_pressure(density):
    density_error = density - target_density
    pressure = density_error * pressure_multiplier
    return pressure

# The average of the pressures of the two particles
def calculate_shared_pressure(densityA, densityB):
    pressureA = convert_density_to_pressure(densityA)
    pressureB = convert_density_to_pressure(densityB)
    return (pressureA + pressureB) / 2


# Returns the coordinate of the cell on the grid
def position_to_cell_coord(particle_pos):
    offset_vector = vector(particle_pos.x % influence_radius, 
                            particle_pos.y % influence_radius, 
                            particle_pos.z % influence_radius)
    cell_coord = (particle_pos - offset_vector) / influence_radius
    return (cell_coord.x, cell_coord.y, cell_coord.z)

def vector_add(p1, p2):
    return (p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2])


# Setting up things before the simulation starts

start_arrangement()

# Update function
def update():
    move_particles()

while True:
    update()
    rate(update_rate)    