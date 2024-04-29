from vpython import *

# Global Variables ---

# Setting up the wall around the Fluid
wall_length = 8.0
wall_width = 4.0
wall_height = 5.0
wall_thickness = 0.1

# Ball parameters
gravity = 9.8
update_rate = 200
delta_time = 0.01

collision_damping = 0.99
ball_radius = 0.4
total_balls = 10


# Instructions to move the scene.
scene.caption = """Right button drag or Ctrl-drag to rotate "camera" to view scene.
To zoom, drag with middle button or Alt/Option depressed, or use scroll wheel.
    On a two-button mouse, middle is left + right.
Shift-drag to pan left/right and up/down.
Touch screen: pinch/extend to zoom, swipe or two-finger rotate."""

# Setting up the scene
wallR = box(pos=vector( wall_length, 0, 0), size=vector(wall_thickness, 2 * wall_height, 2 * wall_width),  color = color.white)
wallL = box(pos=vector(-wall_length, 0, 0), size=vector(wall_thickness, 2 * wall_height, 2 * wall_width),  color = color.white)
wallB = box(pos=vector(0, -wall_height, 0), size=vector(2 * wall_length, wall_thickness, 2 * wall_width),  color = color.white)
wallT = box(pos=vector(0,  wall_height, 0), size=vector(2 * wall_length, wall_thickness, 2 * wall_width),  color = color.white)
wallBK = box(pos=vector(0, 0, -wall_width), size=vector(2 * wall_length, 2 * wall_height, wall_thickness), color = color.white)

# Keeping the balls away from the wall
wall_length = wall_length - wall_thickness * 0.5 - ball_radius
wall_width = wall_width - wall_thickness * 0.5 - ball_radius
wall_height = wall_height - wall_thickness * 0.5 - ball_radius

# Creating the list of balls
balls_list = list()
for i in range(total_balls):
    ball = sphere (color = color.cyan, radius = ball_radius, make_trail=False, retain=200)
    ball.mass = 1.0
    ball.p = vector (-0.1, -0.23, +0.27)
    balls_list.append(ball)


# Moving the balls
def move_balls(balls):
    for ball in balls:
        add_gravity(ball)
        resolve_collisions(ball, collision_damping)


# Adding gravity to the ball's movement
def add_gravity(ball):
    ball.p.y += -gravity * delta_time
    ball.pos.y += ball.p.y * delta_time


# Reflects the ball off if the ball collides the box
def resolve_collisions(ball, collision_damping):
    if not (wall_length > ball.pos.x > -wall_length):
        ball.p.x *= -1 * collision_damping
    if not (wall_height > ball.pos.y > -wall_height):
        ball.p.y *= -1 * collision_damping
    if not (wall_width > ball.pos.z > -wall_width):
        ball.p.z *= -1 * collision_damping
    
    
# Update function
def update():
    move_balls(ball)

while True:
    update()
    rate(update_rate)    

#dt = 0.3
#def move():
#    rate(200, move)
#    ball.pos = ball.pos + (ball.p/ball.mass)*dt
#    if not (side > ball.pos.x > -side):
#        ball.p.x = -ball.p.x
#    if not (side > ball.pos.y > -side):
#        ball.p.y = -ball.p.y
#    if not (side > ball.pos.z > -side):
#        ball.p.z = -ball.p.z
#
#move()