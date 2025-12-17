#!/usr/bin/env python3
"""
Calculate position angle towards celestial equator from CD matrix values
For use with REDUC and similar software that use equatorial reference

This script properly calculates position angles
Reads CD matrix directly from FITS file headers
"""

import math
import sys
import os

try:
    from astropy.io import fits
    ASTROPY_AVAILABLE = True
except ImportError:
    ASTROPY_AVAILABLE = False

def clean_filename(filename):
    """
    Clean filename by removing surrounding quotes and whitespace
    
    Args:
        filename: Raw filename input
    
    Returns:
        Cleaned filename string
    """
    # Strip whitespace first
    filename = filename.strip()
    
    # Remove surrounding double quotes if present
    if filename.startswith('"') and filename.endswith('"'):
        filename = filename[1:-1]
    
    # Remove surrounding single quotes if present
    if filename.startswith("'") and filename.endswith("'"):
        filename = filename[1:-1]
    
    return filename

def read_cd_matrix_from_fits(filename):
    """
    Read CD matrix values from FITS file header
    
    Args:
        filename: Path to FITS file
    
    Returns:
        tuple: (cd1_1, cd1_2, cd2_1, cd2_2) or None if not found
    """
    if not ASTROPY_AVAILABLE:
        print("Error: astropy library is required to read FITS files.")
        print("Install with: pip install astropy")
        return None
    
    try:
        with fits.open(filename) as hdul:
            header = hdul[0].header
            
            # Try to read CD matrix values
            cd1_1 = header.get('CD1_1')
            cd1_2 = header.get('CD1_2') 
            cd2_1 = header.get('CD2_1')
            cd2_2 = header.get('CD2_2')
            
            if all(val is not None for val in [cd1_1, cd1_2, cd2_1, cd2_2]):
                return cd1_1, cd1_2, cd2_1, cd2_2
            
            # If CD matrix not found, try CDELT/CROTA approach
            cdelt1 = header.get('CDELT1')
            cdelt2 = header.get('CDELT2') 
            crota2 = header.get('CROTA2', 0.0)  # Default to 0 if not present
            
            if cdelt1 is not None and cdelt2 is not None:
                print("CD matrix not found, converting from CDELT/CROTA values...")
                crota2_rad = math.radians(crota2)
                cos_rot = math.cos(crota2_rad)
                sin_rot = math.sin(crota2_rad)
                
                cd1_1 = cdelt1 * cos_rot
                cd1_2 = -cdelt2 * sin_rot
                cd2_1 = cdelt1 * sin_rot
                cd2_2 = cdelt2 * cos_rot
                
                return cd1_1, cd1_2, cd2_1, cd2_2
            
            print("Error: Neither CD matrix nor CDELT/CROTA values found in FITS header")
            return None
            
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found")
        return None
    except Exception as e:
        print(f"Error reading FITS file: {e}")
        return None

def calculate_position_angle_north(cd1_1, cd1_2, cd2_1, cd2_2):
    """
    Calculate position angle towards celestial north from CD matrix
    This uses the standard CROTA2 calculation from FITS WCS
    
    Args:
        cd1_1, cd1_2, cd2_1, cd2_2: CD matrix elements from FITS header
    
    Returns:
        Position angle towards north in degrees (0-360)
    """
    # Standard CROTA2 calculation - angle of declination axis relative to pixel Y-axis
    pa_north_rad = math.atan2(-cd1_2, cd1_1)
    pa_north_deg = math.degrees(pa_north_rad)
    
    # Ensure positive angle (0-360 degrees)
    if pa_north_deg < 0:
        pa_north_deg += 360
    
    return pa_north_deg

def calculate_position_angle_equator(cd1_1, cd1_2, cd2_1, cd2_2):
    """
    Calculate position angle towards celestial equator from CD matrix
    This is the angle from north to the direction of the celestial equator
    
    Args:
        cd1_1, cd1_2, cd2_1, cd2_2: CD matrix elements from FITS header
    
    Returns:
        Position angle towards equator in degrees (0-360)
    """
    # First calculate the position angle towards north
    pa_north = calculate_position_angle_north(cd1_1, cd1_2, cd2_1, cd2_2)
    
    # The celestial equator is perpendicular to the north direction
    # So we add 90 degrees to get the equator direction
    pa_equator = (pa_north + 90) % 360
    
    return pa_equator

def calculate_crota2_legacy(cd1_1, cd1_2):
    """
    Calculate CROTA2 using the legacy method for comparison
    
    Args:
        cd1_1, cd1_2: CD matrix elements
    
    Returns:
        CROTA2 in degrees
    """
    crota2_rad = math.atan2(-cd1_2, cd1_1)
    crota2_deg = math.degrees(crota2_rad)
    
    if crota2_deg < 0:
        crota2_deg += 360
        
    return crota2_deg

def main():
    """Main function to get FITS file input and calculate position angles"""
    print("Position Angle Calculator from CD Matrix")
    print("=" * 50)
    
    # Check if filename provided as command line argument
    if len(sys.argv) > 1:
        filename = clean_filename(sys.argv[1])
    else:
        # Get filename from user input
        raw_filename = input("Enter FITS filename (or 'manual' for manual CD matrix input): ").strip()
        filename = clean_filename(raw_filename)
    
    if filename.lower() == 'manual':
        # Manual input mode (original functionality)
        try:
            print("Enter the CD matrix values from your FITS header:")
            cd1_1 = float(input("CD1_1: "))
            cd1_2 = float(input("CD1_2: "))
            cd2_1 = float(input("CD2_1: "))
            cd2_2 = float(input("CD2_2: "))
        except ValueError:
            print("Error: Please enter valid numeric values for the CD matrix elements.")
            return
    else:
        # Read from FITS file
        if not os.path.exists(filename):
            print(f"Error: File '{filename}' not found")
            return
            
        cd_values = read_cd_matrix_from_fits(filename)
        if cd_values is None:
            return
        
        cd1_1, cd1_2, cd2_1, cd2_2 = cd_values
        print(f"Successfully read CD matrix from: {filename}")
    
    try:
        # Calculate angles
        pa_north = calculate_position_angle_north(cd1_1, cd1_2, cd2_1, cd2_2)
        pa_equator = calculate_position_angle_equator(cd1_1, cd1_2, cd2_1, cd2_2)
        crota2 = calculate_crota2_legacy(cd1_1, cd1_2)
        
        print("\n" + "=" * 50)
        print("RESULTS:")
        print("=" * 50)
        print(f"CD Matrix values:")
        print(f"CD1_1 = {cd1_1}")
        print(f"CD1_2 = {cd1_2}")
        print(f"CD2_1 = {cd2_1}")
        print(f"CD2_2 = {cd2_2}")
        print()
        print(f"Position Angle towards North (CROTA2): {pa_north:.2f}°")
        print(f"Position Angle towards Celestial Equator: {pa_equator:.2f}°")
        print(f"  (For REDUC and similar software)")
        print(f"CROTA2 (legacy calculation): {crota2:.2f}°")
        print()
        print("Note: Astrometry.net's orientation from JPEGs can be unreliable")
        print("      due to loss of FITS header information. Use CD matrix values.")
        
        # Calculate pixel scales for verification
        pixel_scale_x = math.sqrt(cd1_1**2 + cd2_1**2) * 3600  # arcsec/pixel
        pixel_scale_y = math.sqrt(cd1_2**2 + cd2_2**2) * 3600  # arcsec/pixel
        print(f"\nPixel scales:")
        print(f"X-axis: {pixel_scale_x:.3f} arcsec/pixel")
        print(f"Y-axis: {pixel_scale_y:.3f} arcsec/pixel")
        
        # Check for image flip
        determinant = cd1_1 * cd2_2 - cd1_2 * cd2_1
        if determinant < 0:
            print("Note: Image appears to be flipped (negative determinant)")
        
    except KeyboardInterrupt:
        print("\nProgram interrupted by user.")

if __name__ == "__main__":
    main()
