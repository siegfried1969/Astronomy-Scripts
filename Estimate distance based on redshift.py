from astropy.cosmology import FlatLambdaCDM

# Define the cosmology model
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def calculate_cosmology():
    try:
        # Prompt user for input
        z_input = input("Enter the redshift (z) of the object: ")
        z = float(z_input)
        
        if z < 0:
            print("Redshift must be a positive number.")
            return

        # Perform calculations
        # Comoving Distance: It tracks where the quasar resides in the universe right now
        comoving = cosmo.comoving_distance(z)
        # Luminosity Distance: This is the optical fading. It is the "effective" distance the light implies based on how dim it looks.
        luminosity = cosmo.luminosity_distance(z)
        # Lookback Distance: This is the photon's life span. It is the duration the light spent traveling through the cosmic dark
        lookback = cosmo.lookback_time(z)
        
        # Display results
        print(f"\n--- Results for z = {z} ---")
        print(f"Comoving Distance:    {comoving.value:.2f} Mpc")
        print(f"Luminosity Distance:  {luminosity.value:.2f} Mpc")
        print(f"Lookback Time:        {lookback.value:.2f} Gyr")
        
    except ValueError:
        print("Error: Please enter a valid numerical value for redshift.")

if __name__ == "__main__":
    calculate_cosmology()