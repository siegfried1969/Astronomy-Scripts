import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Physical constants and parameters
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M1 = 1.989e30    # Mass of star 1 (kg) - approximately solar mass
M2 = 0.5 * M1    # Mass of star 2 (kg) - half the mass of star 1
a = 1.496e11     # Semi-major axis of orbit (m) - approximately 1 AU

# Calculate center of mass position
# Center of mass from star 1: r_cm = M2 * a / (M1 + M2)
M_total = M1 + M2
r1 = M2 * a / M_total  # Distance from star 1 to center of mass
r2 = M1 * a / M_total  # Distance from star 2 to center of mass

print(f"Star 1 mass: {M1/1.989e30:.1f} solar masses")
print(f"Star 2 mass: {M2/1.989e30:.1f} solar masses")
print(f"Star 1 orbital radius: {r1/1.496e11:.3f} AU")
print(f"Star 2 orbital radius: {r2/1.496e11:.3f} AU")

# Calculate orbital velocities
# For circular orbits: v = sqrt(G * M_total / a) * (other_star_distance / a)
# This gives the velocity needed for circular orbit around the center of mass
v1 = np.sqrt(G * M_total * r1 / a**2)  # Velocity of star 1
v2 = np.sqrt(G * M_total * r2 / a**2)  # Velocity of star 2

# Alternative calculation using reduced mass approach:
# v = sqrt(G * M_total / a) for the relative motion
v_rel = np.sqrt(G * M_total / a)
# Then scale by distance from center of mass
v1_alt = v_rel * r1 / a
v2_alt = v_rel * r2 / a

print(f"Star 1 orbital velocity: {v1/1000:.2f} km/s")
print(f"Star 2 orbital velocity: {v2/1000:.2f} km/s")

# Calculate orbital period using the relative motion
T = 2 * np.pi * a / v_rel
print(f"Orbital period: {T/86400:.2f} days")

# Verify that both stars have the same angular velocity
omega1 = v1 / r1
omega2 = v2 / r2
omega = 2 * np.pi / T
print(f"Angular velocity check - omega1: {omega1:.6e}, omega2: {omega2:.6e}, omega: {omega:.6e}")

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
# Both stars orbit the center of mass with the same angular frequency
for i in range(time_steps):
    # Current angle in orbit
    theta = omega * t[i]
    
    # Star 1 position (starts at +r1 on x-axis from center of mass)
    x1[i] = r1 * np.cos(theta)
    y1[i] = r1 * np.sin(theta)
    
    # Star 2 position (starts at -r2 on x-axis from center of mass)
    # Phase difference of π (180 degrees) keeps stars on opposite sides
    x2[i] = r2 * np.cos(theta + np.pi)
    y2[i] = r2 * np.sin(theta + np.pi)

# Create the plot
fig, ax = plt.subplots(figsize=(12, 10))

# Set up the plot limits and appearance
plot_limit = max(r1, r2) * 1.5
ax.set_xlim(-plot_limit, plot_limit)
ax.set_ylim(-plot_limit, plot_limit)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Distance (m)')
ax.set_title('Binary Star System - Unequal Masses\n' + 
             f'Star 1: {M1/1.989e30:.1f} M☉, Star 2: {M2/1.989e30:.1f} M☉')

# Mark the center of mass
ax.plot(0, 0, 'ko', markersize=8, label='Center of Mass')

# Plot the complete orbits as reference
circle1 = plt.Circle((0, 0), r1, fill=False, linestyle='--', alpha=0.5, color='red')
circle2 = plt.Circle((0, 0), r2, fill=False, linestyle='--', alpha=0.5, color='blue')
ax.add_patch(circle1)
ax.add_patch(circle2)

# Initialize the moving elements
# Make star 1 larger since it's more massive
star1, = ax.plot([], [], 'ro', markersize=15, label=f'Star 1 ({M1/1.989e30:.1f} M☉)')
star2, = ax.plot([], [], 'bo', markersize=10, label=f'Star 2 ({M2/1.989e30:.1f} M☉)')
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
# anim.save('binary_stars_unequal.gif', writer='pillow', fps=20)

# Print some interesting facts about the system
print(f"\nSystem Properties:")
print(f"Total system mass: {M_total/1.989e30:.1f} solar masses")
print(f"Separation distance: {a/1.496e11:.2f} AU")
print(f"Star 1 orbital radius: {r1/1.496e11:.3f} AU")
print(f"Star 2 orbital radius: {r2/1.496e11:.3f} AU")
print(f"Star 1 orbital velocity: {v1/1000:.2f} km/s")
print(f"Star 2 orbital velocity: {v2/1000:.2f} km/s")

# Calculate centripetal accelerations
a1_centripetal = v1**2 / r1
a2_centripetal = v2**2 / r2
print(f"Star 1 centripetal acceleration: {a1_centripetal:.2f} m/s²")
print(f"Star 2 centripetal acceleration: {a2_centripetal:.2f} m/s²")

# Verify our physics: centripetal force should equal gravitational force
# For star 1: F_c1 = M1 * a1_centripetal, F_g1 = G * M1 * M2 / a²
# For star 2: F_c2 = M2 * a2_centripetal, F_g2 = G * M1 * M2 / a²
F_c1 = M1 * a1_centripetal
F_c2 = M2 * a2_centripetal
F_g = G * M1 * M2 / a**2

print(f"\nPhysics check:")
print(f"Star 1 centripetal force: {F_c1:.2e} N")
print(f"Star 2 centripetal force: {F_c2:.2e} N")
print(f"Gravitational force: {F_g:.2e} N")
print(f"Force ratio 1 (should be ~1): {F_c1/F_g:.3f}")
print(f"Force ratio 2 (should be ~1): {F_c2/F_g:.3f}")

# Additional verification: momentum conservation
# Total momentum should be zero (center of mass frame)
p1 = M1 * v1  # Momentum magnitude of star 1
p2 = M2 * v2  # Momentum magnitude of star 2
print(f"\nMomentum check:")
print(f"Star 1 momentum magnitude: {p1:.2e} kg⋅m/s")
print(f"Star 2 momentum magnitude: {p2:.2e} kg⋅m/s")
print(f"Momentum ratio (should be ~1): {p1/p2:.3f}")