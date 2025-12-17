import os
import glob
from astropy.io import fits
import numpy as np


def stack_fits_files(directory):
    # Remove surrounding quotes if present
    directory = directory.strip('"')
    directory = os.path.normpath(directory)

    # Get all FITS files in the directory ending with _00000.fits, _00001.fits, etc.
    fits_files = glob.glob(os.path.join(directory, "*_?????.fits"))
    fits_files.sort()  # Sort files for consistent ordering

    if not fits_files:
        print("Error: No matching FITS files found in the directory.")
        return

    # Read the first file to get the shape and header
    try:
        with fits.open(fits_files[0]) as hdul:
            data_shape = hdul[0].data.shape
            header = hdul[0].header
    except Exception as e:
        print(f"Error reading the first FITS file: {e}")
        return

    # Initialize a 3D array for the cube
    cube = np.empty((len(fits_files), data_shape[0], data_shape[1]), dtype=np.float32)

    # Read and stack all FITS files
    for i, fits_file in enumerate(fits_files):
        try:
            with fits.open(fits_file) as hdul:
                cube[i] = hdul[0].data.astype(np.float32)
        except Exception as e:
            print(f"Error reading {fits_file}: {e}")
            return

    # Generate the output filename by replacing the last 5 characters before .fits with CUBE
    base_name = os.path.basename(fits_files[0])
    output_name = base_name[:-10] + "CUBE.fits"  # Remove _00000 and replace with CUBE
    output_path = os.path.join(directory, output_name)

    # Write the 3D cube to a new FITS file
    try:
        hdu = fits.PrimaryHDU(cube, header=header)
        hdu.writeto(output_path, overwrite=True)
        print(f"Success! Stacked cube saved to: {output_path}")
    except Exception as e:
        print(f"Error writing the output FITS file: {e}")


if __name__ == "__main__":
    directory = input("Enter the directory path containing FITS files: ").strip()
    stack_fits_files(directory)
