import math
from datetime import datetime, timezone
from astroquery.gaia import Gaia


def propagate_coordinates(ra_2016, dec_2016, pmra, pmdec, epoch_2016, target_epoch):
    """
    Propagate RA and DEC from epoch_2016 to target_epoch.

    Args:
        ra_2016 (float): RA at J2016.0 (degrees)
        dec_2016 (float): DEC at J2016.0 (degrees)
        pmra (float): Proper motion in RA (mas/yr)
        pmdec (float): Proper motion in DEC (mas/yr)
        epoch_2016 (float): Reference epoch (2016.0)
        target_epoch (float): Target epoch (e.g., current year as a decimal)

    Returns:
        tuple: (ra_target, dec_target) in degrees
    """
    delta_t = target_epoch - epoch_2016

    # Convert proper motion from mas/yr to deg/yr
    pmra_deg_per_yr = pmra / (3600.0 * 1000.0)
    pmdec_deg_per_yr = pmdec / (3600.0 * 1000.0)

    # Calculate new DEC
    dec_target = dec_2016 + pmdec_deg_per_yr * delta_t

    # Calculate new RA (accounting for cos(DEC))
    cos_dec = math.cos(math.radians(dec_2016))
    ra_target = ra_2016 + (pmra_deg_per_yr * delta_t) / cos_dec

    return ra_target, dec_target


def deg_to_hms(ra_deg):
    """Convert RA from degrees to hours, minutes, seconds."""
    ra_h = ra_deg / 15.0
    h = int(ra_h)
    m = int((ra_h - h) * 60)
    s = ((ra_h - h) * 60 - m) * 60
    return h, m, s


def deg_to_dms(dec_deg):
    """Convert DEC from degrees to degrees, arcminutes, arcseconds."""
    dec_deg_abs = abs(dec_deg)
    d = int(dec_deg_abs)
    m = int((dec_deg_abs - d) * 60)
    s = ((dec_deg_abs - d) * 60 - m) * 60
    if dec_deg < 0:
        d = -d
    return d, m, s


# Prompt the user for a Gaia DR3 source ID
source_id = input("Enter the Gaia DR3 source ID (e.g., 3346233862409398272): ").strip()
# Example: source_id = 3346233862409398272

# Query Gaia DR3 for the source ID
query = f"""
SELECT
    ra, dec, pmra, pmdec, ref_epoch
FROM gaiadr3.gaia_source
WHERE source_id = {source_id}
"""
job = Gaia.launch_job(query)
result = job.get_results()

# Extract data from the query result
ra_2016 = result["ra"][0]
dec_2016 = result["dec"][0]
pmra = result["pmra"][0]
pmdec = result["pmdec"][0]
epoch_2016 = result["ref_epoch"][0]

# Calculate J2000 coordinates
ra_j2000, dec_j2000 = propagate_coordinates(
    ra_2016, dec_2016, pmra, pmdec, epoch_2016, 2000.0
)

# Get current UTC time and convert to decimal year
now = datetime.now(timezone.utc)
current_year = (
    now.year
    + (
        now.timetuple().tm_yday
        + now.hour / 24.0
        + now.minute / 1440.0
        + now.second / 86400.0
    )
    / 366.0
)

# Calculate current coordinates
ra_now, dec_now = propagate_coordinates(
    ra_2016, dec_2016, pmra, pmdec, epoch_2016, current_year
)

# Convert to human-readable format
ra_j2000_h, ra_j2000_m, ra_j2000_s = deg_to_hms(ra_j2000)
dec_j2000_d, dec_j2000_m, dec_j2000_s = deg_to_dms(dec_j2000)

ra_now_h, ra_now_m, ra_now_s = deg_to_hms(ra_now)
dec_now_d, dec_now_m, dec_now_s = deg_to_dms(dec_now)

# Print results
print(f"J2000 RA: {ra_j2000_h}h {ra_j2000_m}m {ra_j2000_s:.2f}s")
print(f"J2000 DEC: {dec_j2000_d}° {dec_j2000_m}' {dec_j2000_s:.2f}\"")
print(f"Current RA (J{current_year:.6f}): {ra_now_h}h {ra_now_m}m {ra_now_s:.2f}s")
print(f"Current DEC (J{current_year:.6f}): {dec_now_d}° {dec_now_m}' {dec_now_s:.2f}\"")
