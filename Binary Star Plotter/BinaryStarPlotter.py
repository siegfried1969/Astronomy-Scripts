import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Physical constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M = 1.989e30     # Mass of each star (kg) - approximately solar mass
a = 1.496e11     # Semi-major axis of orbit (m) - approximately 1 AU

# Calculate orbital parameters for identical binary stars
# For identical masses in circular orbit, each star orbits around their center of mass
# The center of mass is exactly halfway between them
# Each star is at distance a/2 from the center of mass
r = a / 2  # Distance from center of mass to each star

# Calculate orbital velocity using Kepler's laws
# For circular orbit: v = sqrt(GM/r_orbit)
# But in binary system, we need to account for both masses
# The reduced mass system gives us: v = sqrt(G(M1+M2)/a) for relative motion
# Since each star moves in circle of radius r = a/2 around center of mass:
# v = sqrt(GM/r) where the effective mass is M (the other star)
v = np.sqrt(G * M / r)

# Calculate orbital period using T = 2πr/v
T = 2 * np.pi * r / v
print(f"Orbital period: {T/86400:.2f} days")
print(f"Orbital velocity: {v/1000:.2f} km/s")

# Set up the simulation parameters
dt = T / 1000  # Time step (1/1000 of orbital period for smooth animation)
t_max = 2 * T  # Total simulation time (2 complete orbits)
time_steps = int(t_max / dt)

# Initialize arrays to store positions
t = np.linspace(0, t_max, time_steps)
x1 = np.zeros(time_steps)  # x position of star 1
y1 = np.zeros(time_steps)  # y position of star 1
x2 = np.zeros(time_steps)  # x position of star 2
y2 = np.zeros(time_steps)  # y position of star 2

# Calculate positions over time
# Angular frequency: ω = v/r = 2π/T
omega = 2 * np.pi / T

for i in range(time_steps):
    # Current angle in orbit
    theta = omega * t[i]
    
    # Star 1 position (starts at +r on x-axis)
    x1[i] = r * np.cos(theta)
    y1[i] = r * np.sin(theta)
    
    # Star 2 position (starts at -r on x-axis, opposite side)
    # Phase difference of π (180 degrees) keeps stars on opposite sides
    x2[i] = r * np.cos(theta + np.pi)
    y2[i] = r * np.sin(theta + np.pi)

# Create the plot
fig, ax = plt.subplots(figsize=(10, 10))

# Set up the plot limits and appearance
plot_limit = r * 1.5
ax.set_xlim(-plot_limit, plot_limit)
ax.set_ylim(-plot_limit, plot_limit)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_title('Binary Star System - Circular Orbit\n(Two identical stars orbiting their center of mass)')

# Mark the center of mass
ax.plot(0, 0, 'ko', markersize=8, label='Center of Mass')

# Plot the complete orbits as reference
circle1 = plt.Circle((0, 0), r, fill=False, linestyle='--', alpha=0.5, color='blue')
ax.add_patch(circle1)

# Initialize the moving elements
star1, = ax.plot([], [], 'ro', markersize=12, label='Star 1')
star2, = ax.plot([], [], 'bo', markersize=12, label='Star 2')
trail1, = ax.plot([], [], 'r-', alpha=0.3, linewidth=1)
trail2, = ax.plot([], [], 'b-', alpha=0.3, linewidth=1)

# Add connecting line to show the system
connection, = ax.plot([], [], 'k--', alpha=0.5, linewidth=1)

ax.legend()

# Animation function
def animate(frame):
    """
    Animation function called for each frame
    Updates positions of stars and their trails
    """
    # Length of trail to show (last 100 points)
    trail_length = min(100, frame + 1)
    start_idx = max(0, frame - trail_length + 1)
    
    # Update star positions
    star1.set_data([x1[frame]], [y1[frame]])
    star2.set_data([x2[frame]], [y2[frame]])
    
    # Update trails
    trail1.set_data(x1[start_idx:frame+1], y1[start_idx:frame+1])
    trail2.set_data(x2[start_idx:frame+1], y2[start_idx:frame+1])
    
    # Update connection line between stars
    connection.set_data([x1[frame], x2[frame]], [y1[frame], y2[frame]])
    
    return star1, star2, trail1, trail2, connection

# Create and run animation
# interval controls animation speed (milliseconds between frames)
anim = animation.FuncAnimation(fig, animate, frames=time_steps, 
                             interval=50, blit=True, repeat=True)

# Display the plot
plt.tight_layout()
plt.show()

# Optional: Save as GIF (uncomment the line below)
# anim.save('binary_stars.gif', writer='pillow', fps=20)

# Print some interesting facts about the system
print(f"\nSystem Properties:")
print(f"Total system mass: {2*M/1.989e30:.1f} solar masses")
print(f"Separation distance: {a/1.496e11:.2f} AU")
print(f"Each star's orbital radius: {r/1.496e11:.2f} AU")
print(f"Orbital velocity: {v/1000:.2f} km/s")
print(f"Centripetal acceleration: {v**2/r:.2f} m/s²")

# Verify our physics: centripetal force should equal gravitational force
# Centripetal force: F_c = M * v² / r
# Gravitational force: F_g = G * M² / (2r)² = G * M² / (4r²)
# For circular orbit, these should be equal
F_centripetal = M * v**2 / r
F_gravitational = G * M**2 / (2*r)**2
print(f"\nPhysics check:")
print(f"Centripetal force: {F_centripetal:.2e} N")
print(f"Gravitational force: {F_gravitational:.2e} N")
print(f"Ratio (should be ~1): {F_centripetal/F_gravitational:.3f}")
