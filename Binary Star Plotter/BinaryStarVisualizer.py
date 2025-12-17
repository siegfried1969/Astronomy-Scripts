import matplotlib.pyplot as plt
import numpy as np

def visualize_double_star(data_string):
    """
    Parses double star observation data and visualizes the secondary star's position
    relative to the primary star (at the origin).

    The input string must contain 9 space-separated values (Year, PA, SEP) for 
    three observations. The Separation (SEP) values MUST be in arcseconds.
    
    Args:
        data_string (str): Input string, e.g., "1834 302.00 3.00 2014 298.00 3.37 2025 297.07 3.40"
    """
    
    # 1. Parsing the Input String
    tokens = data_string.split()
    
    if len(tokens) != 9:
        print("Error: Input string must contain exactly 9 space-separated values (Year, PA, SEP for 3 observations).")
        return

    observations = []
    for i in range(0, 9, 3):
        try:
            year = int(tokens[i])
            pa_deg = float(tokens[i+1])
            # Renamed variable here for clarity, though it was already correct:
            sep_arcsec = float(tokens[i+2]) 
            observations.append((year, pa_deg, sep_arcsec))
        except ValueError as e:
            print(f"Error parsing token: {e}. Check that year is integer and PA/SEP are floats.")
            return

    # Data conversion and Cartesian calculation
    years = []
    x_positions = []
    y_positions = []

    # The loop is now simple: data is already in the correct units (degrees, arcseconds).
    for year, pa_deg, sep_arcsec in observations:
        # **NO CONVERSION NEEDED: SEP is already in arcseconds**

        # Convert Position Angle (PA) from degrees to radians
        pa_rad = np.deg2rad(pa_deg)

        # Calculate Cartesian coordinates (x, y)
        # Astronomical convention: PA=0 is North (Y+), PA=90 is East (X+)
        x = sep_arcsec * np.sin(pa_rad)
        y = sep_arcsec * np.cos(pa_rad)

        years.append(year)
        x_positions.append(x)
        y_positions.append(y)
        
    # ... (Rest of the plotting code is correct and unchanged)
    # ... (Plotting code omitted for brevity)
    
    # 2. Visualization Setup
    
    # Set up the figure and axes with a square aspect ratio
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Define the maximum scale as requested (400 arcseconds)
    max_scale = 400

    # Plot the Primary Star at the center (0, 0)
    ax.plot(0, 0, '*', color='gold', markersize=20, markeredgecolor='black', label='Primary Star (0, 0)')

    # Plot Secondary Star positions, using year for color mapping
    scatter = ax.scatter(
        x_positions, 
        y_positions, 
        s=30, 
        color='black',
        zorder=5,
        edgecolors='white'
    )

    # 3. Grid and Aesthetics
    
    # Set axis limits to 400 arcseconds. 
    # Reverse X-axis (West to Right, East to Left) for standard astronomical charting.
    ax.set_xlim(max_scale, -max_scale) 
    ax.set_ylim(-max_scale, max_scale) # North is up, South is down

    # Ensure the plot maintains a 1:1 aspect ratio
    ax.set_aspect('equal', adjustable='box')
    
    # Add concentric separation rings for the "visual grid"
    # Rings now go up to 400 arcseconds
    for r in [100, 200, 300, 400]:
        # Draw a thin, dashed circle
        circle = plt.Circle((0, 0), r, color='gray', fill=False, linestyle='--', linewidth=0.7, zorder=0, alpha=0.5)
        ax.add_artist(circle)
        # Label the ring separation on the right side
        ax.text(r, 0, f'{r}"', color='gray', ha='left', va='center', fontsize=8)

    # Add primary axes lines
    ax.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax.axvline(0, color='black', linestyle='-', linewidth=0.5)
    
    # Label North, South, East, West directions
    offset = 30
    ax.text(0, max_scale + offset, 'NORTH (0°)', ha='center', va='center', fontsize=12, color='navy', weight='bold')
    ax.text(0, -max_scale - offset, 'SOUTH (180°)', ha='center', va='center', fontsize=12, color='navy', weight='bold')
    ax.text(max_scale + offset, 0, 'EAST (90°)', ha='center', va='center', fontsize=12, color='navy', weight='bold')
    ax.text(-max_scale - offset, 0, 'WEST (270°)', ha='center', va='center', fontsize=12, color='navy', weight='bold')
    
    # Set final labels and title
    offset = 30
    ax.set_title('Double Star Observation Plot', fontsize=16, weight='bold', pad=offset)

    # Remove default tick labels on the axes, as we use N/S/E/W labels
    ax.set_xticks([])
    ax.set_yticks([])
    # Clear the previously set title and place it explicitly above the 'NORTH' label
    ax.set_title('') # remove the default title (we'll place a data-coordinate title)

    # Compute title y-coordinate so it sits above the NORTH text
    # 'offset' was previously set to 30 when the direction labels were placed.
    title_gap = 12 # additional gap above the NORTH label in data units
    title_y = max_scale + offset + title_gap

    ax.text(
        0, title_y,
        'Double Star Observation Plot with Path',
        ha='center',
        va='bottom',
        fontsize=16,
        weight='bold',
        color='black',
        zorder=10,
        clip_on=False
    )   
    plt.grid(False) # Turn off Matplotlib's default grid, as we drew our own circles/lines
    plt.tight_layout()
    
    print("--- Plotting function successfully executed. Attempting to display visualization. ---")
    plt.show()

if __name__ == '__main__':
    # Example data following the specified format:
    # Year1 PA1 SEP1 Year2 PA2 SEP2 Year3 PA3 SEP3
    
    # Note: SEP values MUST be in arcseconds ("). The plot limit is 400 arcseconds.
    example_input = input()
    # E.g., for stars at ~140 arcseconds separation: "1895 352.00 123.50 2014 358.00 140.70 2025 358.90 142.37" 
    
    
    print("Running double star visualizer with example data:")
    visualize_double_star(example_input)