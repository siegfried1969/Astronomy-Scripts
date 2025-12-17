import sys
import math
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.utils.exceptions import AstropyWarning
import warnings

# --- WARNING SUPPRESSION ---
warnings.filterwarnings("ignore", category=AstropyWarning)
warnings.filterwarnings("ignore", category=u.UnitsWarning)

# --- CONFIGURATION ---
TARGET_TABLE = "gaiadr3.gaia_source"
COLUMNS = (
    "source_id, ra, dec, parallax, parallax_error, "
    "pmra, pmdec, pmra_error, pmdec_error, "
    "phot_g_mean_mag, ruwe"
)


def get_user_input():
    print("\n--- Gaia DR3 Binary Analyzer (Law of Cosines) ---")
    # Default Example: 61 Cygni
    # ID 1: 1872066417637336192
    # ID 2: 1872066421962376192

    p_id = input("Enter PRIMARY ID:   ").strip()
    s_id = input("Enter SECONDARY ID: ").strip()

    if not p_id or not s_id:
        sys.exit("Error: Both IDs required.")
    if p_id == s_id:
        sys.exit("Error: Identical IDs entered.")

    return int(p_id), int(s_id)


def fetch_pair_data(id1, id2):
    print(f"\n... Fetching raw data from Gaia DR3 ...")
    query = f"SELECT {COLUMNS} FROM {TARGET_TABLE} WHERE source_id IN ({id1}, {id2})"

    try:
        job = Gaia.launch_job(query)
        results = job.get_results()

        if len(results) < 2:
            return None, None

        def clean_row(row):
            data = {}
            for col in results.colnames:
                val = row[col]
                if np.ma.is_masked(val):
                    data[col] = np.nan
                elif isinstance(val, u.Quantity):
                    data[col] = val.value
                else:
                    data[col] = val
            return data

        rows_by_id = {row["source_id"]: row for row in results}
        return clean_row(rows_by_id[id1]), clean_row(rows_by_id[id2])

    except Exception as e:
        print(f"Query Error: {e}")
        return None, None


def print_raw_data(label, s):
    """Prints the raw inputs from Gaia before any math is done."""
    print(f"\n--- {label} Raw Data ---")
    print(f"  ID:       {s['source_id']}")
    print(f"  RA/Dec:   {s['ra']:.5f}, {s['dec']:.5f} deg")
    print(
        f"  Parallax: {s['parallax']:.5f} \u00b1 {s['parallax_error']:.5f} mas"
    )  # Using \u00B1 for ±
    print(f"  PM:       RA {s['pmra']:.2f}, Dec {s['pmdec']:.2f} mas/yr")
    print(f"  RUWE:     {s['ruwe']:.3f} (Fit Quality)")


def calculate_physics_verbose(s1, s2):
    print(f"\n{'='*20} STEP-BY-STEP CALCULATION {'='*20}")

    # --- STEP 1: INVERSE PARALLAX (Distance) ---
    if np.isnan(s1["parallax"]) or s1["parallax"] <= 0:
        print("Error: Invalid parallax. Cannot calculate distance.")
        return

    # Math: d = 1000 / parallax
    d1_pc = 1000.0 / s1["parallax"]
    d2_pc = 1000.0 / s2["parallax"]

    # Math: error = d * (p_error / p)
    d1_err = d1_pc * (s1["parallax_error"] / s1["parallax"])
    d2_err = d2_pc * (s2["parallax_error"] / s2["parallax"])

    print("\n[Step 1: Distance Calculation]")
    print(f"  Formula: Dist = 1000 / Parallax_mas")
    print(
        f"  Star 1:  1000 / {s1['parallax']:.4f} = {d1_pc:.4f} pc (\u00b1{d1_err:.4f})"
    )
    print(
        f"  Star 2:  1000 / {s2['parallax']:.4f} = {d2_pc:.4f} pc (\u00b1{d2_err:.4f})"
    )

    # --- STEP 2: TRUE 3D SPATIAL SEPARATION (Law of Cosines) ---
    # 1. Calculate Angular Separation
    c1 = SkyCoord(ra=s1["ra"] * u.deg, dec=s1["dec"] * u.deg, frame="icrs")
    c2 = SkyCoord(ra=s2["ra"] * u.deg, dec=s2["dec"] * u.deg, frame="icrs")
    sep_arcsec = c1.separation(c2).arcsec

    # 2. Convert Distances to AU
    AU_PER_PC = 206265.0
    a_au = d1_pc * AU_PER_PC
    b_au = d2_pc * AU_PER_PC

    # 3. Convert Angle to Radians
    theta_rad = math.radians(sep_arcsec / 3600.0)

    # 4. Law of Cosines: c^2 = a^2 + b^2 - 2ab cos(theta)
    c_squared = (a_au**2) + (b_au**2) - (2 * a_au * b_au * math.cos(theta_rad))
    spatial_sep_au = math.sqrt(c_squared)

    print("\n[Step 2: True 3D Separation (Law of Cosines)]")
    print(f"  Angular Separation: {sep_arcsec:.4f} arcsec")
    print(f"  Theta (Radians):    {theta_rad:.9f} rad")
    print(f"  Distance 1 (a):     {a_au:,.0f} AU")
    print(f"  Distance 2 (b):     {b_au:,.0f} AU")
    print(f"  Formula: c = sqrt(a^2 + b^2 - 2ab cos(theta))")
    print(f"  Calculated 3D Separation: {spatial_sep_au:,.2f} AU")

    # --- STEP 3: NUGENT CRITERION A (Overlap) ---
    delta_dist_pc = abs(d1_pc - d2_pc)
    combined_error_pc = d1_err + d2_err
    is_overlapping = delta_dist_pc <= combined_error_pc

    # Added flush=True to prevent repeated terminal output
    print("\n[Step 3: The Overlap Check (Criterion A)]", flush=True)
    print(
        "  Logic: Is the gap between stars smaller than their error bars combined?",
        flush=True,
    )
    print(f"  Gap (|d1 - d2|):    {delta_dist_pc:.4f} pc", flush=True)
    print(f"  Max Allowed (e1+e2): {combined_error_pc:.4f} pc", flush=True)
    print(
        f"  Result: {delta_dist_pc:.4f} <= {combined_error_pc:.4f} ? -> {is_overlapping}",
        flush=True,
    )
    if is_overlapping:
        print("  => Distances OVERLAP within errors.", flush=True)
        # --- STEP 4: DETERMINE MINIMAL SEPARATION USING THE OVERLAP CRITERION---
        # 1. Determine the closest common distance assuming overlap
        d1_min_dist = d1_pc - d1_err
        d2_min_dist = d2_pc - d2_err
        min_comm_distance_pc = max(d1_min_dist, d2_min_dist)
        min_possible_sep_au = sep_arcsec * min_comm_distance_pc  # Small angle approx
        print("\n[Step 4: Minimum Possible Separation calculation]", flush=True)
        print(
            f"  Mimimum possible separation using the overlap criterium: ({min_possible_sep_au:.0f} AU)",
            flush=True,
        )

    # --- FINAL DECISION ---
    print(f"\n{'='*20} FINAL VERDICT {'='*20}", flush=True)

    # Case 1: The stars are statistically at different distances
    if not is_overlapping:
        print(">> OPTICAL DOUBLE (Chance Alignment)", flush=True)
        print(
            f"  Why: The stars are {delta_dist_pc:.2f} pc apart in depth.", flush=True
        )
        print(
            f"  The 3D separation is {spatial_sep_au:,.0f} AU, which is physically impossible for a binary.",
            flush=True,
        )
        return

    # Case 2: The stars overlap in distance (could be a binary)
    if min_possible_sep_au <= 2000:
        print(">> PHYSICAL BINARY (Likely)", flush=True)
        print(
            f"  Why: Distance errors overlap AND projected separation ({min_possible_sep_au:.0f} AU) is < 2,000 AU.",
            flush=True,
        )

        # PM Check
        pm_diff = np.sqrt(
            (s1["pmra"] - s2["pmra"]) ** 2 + (s1["pmdec"] - s2["pmdec"]) ** 2
        )
        print(f"  PM Difference: {pm_diff:.2f} mas/yr", flush=True)
    else:
        print(">> COMMON PROPER MOTION / WIDE BINARY", flush=True)
        print(
            f"  Why: Distances overlap, but separation ({min_possible_sep_au:.0f} AU) > 2,000 AU.",
            flush=True,
        )


if __name__ == "__main__":
    id1, id2 = get_user_input()
    data1, data2 = fetch_pair_data(id1, id2)

    if data1 and data2:
        print_raw_data("Primary", data1)
        print_raw_data("Secondary", data2)
        calculate_physics_verbose(data1, data2)
    else:
        print("\nError: Could not retrieve data.")
