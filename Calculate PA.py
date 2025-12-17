import numpy as np


def hms_to_deg(ra_h, ra_m, ra_s):
    """Converts RA from H:M:S format to Decimal Degrees."""
    return 15.0 * (ra_h + ra_m / 60.0 + ra_s / 3600.0)


def dms_to_deg(dec_d, dec_m, dec_s, sign):
    """Converts Dec from D:M:S format to Decimal Degrees."""
    dec_deg = dec_d + dec_m / 60.0 + dec_s / 3600.0
    return dec_deg * sign


def get_star_input(star_name):
    """Prompts user for RA, Dec, and Pixel coordinates for a star."""
    print(f"\n--- Input for {star_name} ---")

    # RA Input
    ra_h = float(input(f"Enter {star_name} Right Ascension (Hours): "))
    ra_m = float(input(f"Enter {star_name} Right Ascension (Minutes): "))
    ra_s = float(input(f"Enter {star_name} Right Ascension (Seconds): "))

    # Dec Input
    dec_d = float(input(f"Enter {star_name} Declination (Degrees, ABS value): "))
    dec_m = float(input(f"Enter {star_name} Declination (Minutes): "))
    dec_s = float(input(f"Enter {star_name} Declination (Seconds): "))
    dec_sign_str = input(f"Enter {star_name} Declination Sign (+ or -): ")
    dec_sign = 1 if dec_sign_str.strip() in ('+', '1') else -1

    # Pixel Input
    x_coord = float(input(f"Enter {star_name} FITS X-coordinate (pixels): "))
    y_coord = float(input(f"Enter {star_name} FITS Y-coordinate (pixels): "))

    return ra_h, ra_m, ra_s, dec_d, dec_m, dec_s, dec_sign, x_coord, y_coord


def calculate_camera_orientation():
    """
    Guides the user through input and calculates the Camera Angle (Phi).
    """
    print("\n--- Camera Orientation Calculator for Astrometry ---")

    # Get all inputs
    ra1_h, ra1_m, ra1_s, dec1_d, dec1_m, dec1_s, dec1_sign, x1, y1 = get_star_input("Star A (Reference)")
    ra2_h, ra2_m, ra2_s, dec2_d, dec2_m, dec2_s, dec2_sign, x2, y2 = get_star_input("Star B (Companion)")

    pixel_scale = float(input("\nEnter Pixel Scale (arcseconds/pixel): "))

    # --- 1. Converting Coordinates ---
    print("\n--- 1. Converting Coordinates ---")

    # Convert RA/Dec to Decimal Degrees
    alpha1 = hms_to_deg(ra1_h, ra1_m, ra1_s)
    delta1 = dms_to_deg(dec1_d, dec1_m, dec1_s, dec1_sign)
    alpha2 = hms_to_deg(ra2_h, ra2_m, ra2_s)
    delta2 = dms_to_deg(dec2_d, dec2_m, dec2_s, dec2_sign)

    # Convert Decimal Degrees to Radians for trigonometric functions
    delta1_rad = np.deg2rad(delta1)

    print(f"Star A: RA={alpha1:.6f} deg, Dec={delta1:.6f} deg")
    print(f"Star B: RA={alpha2:.6f} deg, Dec={delta2:.6f} deg")

    # --- 2. Calculate True Celestial Position Angle (PA_cel) ---
    print("\n--- 2. Calculating Celestial PA (PA_cel) ---")

    delta_delta = delta2 - delta1

    # Scaled Delta Right Ascension (planar approximation)
    delta_alpha_prime = (alpha2 - alpha1) * np.cos(delta1_rad)

    # PA_cel: East from North
    pa_cel_rad = np.arctan2(delta_alpha_prime, delta_delta)
    pa_cel = np.rad2deg(pa_cel_rad)

    # Ensure PA is in the 0 to 360 range
    if pa_cel < 0:
        pa_cel += 360.0

    print(f"PA_cel (True Celestial Angle, B wrt A): {pa_cel:.2f} degrees")

    # --- 3. Calculate Pixel Position Angle (theta_pix) ---
    print("\n--- 3. Calculating Pixel PA (theta_pix) ---")

    delta_x = x2 - x1
    delta_y = y2 - y1

    # theta_pix: Angle relative to the positive Y-axis of the image.
    theta_pix_rad = np.arctan2(delta_x, delta_y)
    theta_pix = np.rad2deg(theta_pix_rad)

    # Ensure theta_pix is in the 0 to 360 range
    if theta_pix < 0:
        theta_pix += 360.0

    print(f"Delta X (pixels): {delta_x:.2f}")
    print(f"Delta Y (pixels): {delta_y:.2f}")
    print(f"theta_pix (Pixel Angle, B wrt A): {theta_pix:.2f} degrees")

    # --- 4. Calculate Camera Angle (Phi) ---
    print("\n--- 4. Calculating Camera Angle (Phi) ---")

    # Phi is the angle required to rotate the pixel system to match the celestial system
    camera_angle_phi = pa_cel - theta_pix

    # Normalize Phi to the 0 to 360 range
    camera_angle_phi = camera_angle_phi % 360.0

    print(f"FINAL RESULT: Camera Angle (Phi) = {camera_angle_phi:.2f} degrees")

    print("\n--- Interpretation (Y-axis direction) ---")
    if camera_angle_phi > 180 and camera_angle_phi <= 270:
        rotation_from_south = camera_angle_phi - 180.0
        print(f"The camera's positive Y-axis points {rotation_from_south:.2f} degrees East of South.")
    elif camera_angle_phi > 270:
        rotation_from_north = 360.0 - camera_angle_phi
        print(f"The camera's positive Y-axis points {rotation_from_north:.2f} degrees West of North.")
    elif camera_angle_phi <= 180 and camera_angle_phi > 90:
        rotation_from_south = 180.0 - camera_angle_phi
        print(f"The camera's positive Y-axis points {rotation_from_south:.2f} degrees West of South.")
    else:  # 0 to 90
        rotation_from_north = camera_angle_phi
        print(f"The camera's positive Y-axis points {rotation_from_north:.2f} degrees East of North.")

    return camera_angle_phi


if __name__ == '__main__':
    # Run the main function
    calculate_camera_orientation()