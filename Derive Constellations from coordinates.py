"""
Script to map astronomical coordinates to their corresponding constellations.

Requirements:
- Python 3.x
- astropy package (install with: pip install astropy)

Usage:
1. Save this script as coordinate_to_constellation.py
2. Run the script and provide the input file path when prompted
3. The output will be saved to a file named "coordinates_with_constellations.txt" in the same directory
"""

import os
from astropy.coordinates import SkyCoord
from astropy import units as u

def parse_coordinate(coord_str):
    """
    Parse a coordinate string in format HHMMSS.ss+DDMMSS.s or HHMMSS.ss-DDMMSS.s
    into RA and Dec components.

    Returns:
        tuple: (ra_str, dec_str) where ra_str is in HH:MM:SS.ss and dec_str is in ±DD:MM:SS.s
    """
    # Split into RA and Dec parts
    if '+' in coord_str:
        ra_part, dec_part = coord_str.split('+')
        dec_sign = '+'
    else:
        ra_part, dec_part = coord_str.split('-')
        dec_sign = '-'

    # Format RA as HH:MM:SS.ss
    ra_h = ra_part[:2]
    ra_m = ra_part[2:4]
    ra_s = ra_part[4:]
    ra_str = f"{ra_h}:{ra_m}:{ra_s}"

    # Format Dec as ±DD:MM:SS.s
    dec_d = dec_part[:2]
    dec_m = dec_part[2:4]
    dec_s = dec_part[4:]
    dec_str = f"{dec_sign}{dec_d}:{dec_m}:{dec_s}"

    return ra_str, dec_str

def get_constellation(ra_str, dec_str):
    """
    Get the constellation name for given RA and Dec coordinates.

    Args:
        ra_str (str): Right Ascension in HH:MM:SS.ss format
        dec_str (str): Declination in ±DD:MM:SS.s format

    Returns:
        str: Constellation name
    """
    try:
        # Create a SkyCoord object with proper units
        coord = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg), frame='icrs')
        # Get the constellation name
        return coord.get_constellation()
    except Exception as e:
        print(f"Error processing coordinates {ra_str}, {dec_str}: {str(e)}")
        return "Unknown"

def process_file(input_path):
    """
    Process the input file and write results to output file.

    Args:
        input_path (str): Path to input text file
    """
    try:
        # Get the directory of the input file
        input_dir = os.path.dirname(input_path)
        # Create output file path in the same directory
        output_path = os.path.join(input_dir, "coordinates_with_constellations.txt")

        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                try:
                    # Parse the coordinate
                    ra_str, dec_str = parse_coordinate(line)
                    # Get constellation
                    constellation = get_constellation(ra_str, dec_str)
                    # Write to output file
                    outfile.write(f"{line} ; {constellation}\n")
                except Exception as e:
                    print(f"Error processing line '{line}': {str(e)}")
                    outfile.write(f"{line} ; Error\n")

        print(f"\nProcessing complete!")
        print(f"Results saved to: {output_path}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found.")
    except IOError as e:
        print(f"Error reading/writing files: {str(e)}")

if __name__ == "__main__":
    print("Astronomical Coordinate to Constellation Mapper")
    print("----------------------------------------------")

    # Prompt user for input file path
    input_path = input("Enter the path to your coordinates file: ").strip()

    # Remove quotes if present (in case user copies path with quotes)
    input_path = input_path.strip('"\'')

    # Verify the file exists
    if not os.path.isfile(input_path):
        print(f"Error: The file '{input_path}' does not exist.")
        exit(1)

    print(f"\nProcessing file: {input_path}")
    process_file(input_path)