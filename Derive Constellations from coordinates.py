#!/usr/bin/env python3
"""
Append the IAU constellation name to each coordinate in a text file.

Input line format : HHMMSS.ss+DDMMSS.s   (e.g. 092118.79-031225.4)
Output line format: <coordinate_string> ; <constellation_name>

Setup:
    pip install astropy numpy

Usage:
    python add_constellation.py                       # prompts for the source file path
    python add_constellation.py Coordinates.txt       # path given on the command line
    python add_constellation.py Coordinates.txt out.txt --abbr   # 3-letter IAU abbreviations

If no input path is given, the script asks for one and keeps asking until it
points to an existing file. The output defaults to <input>_out.txt, next to
the source file.

Notes:
  * Coordinates are assumed to be J2000/ICRS. astropy converts them internally to
    the B1875 frame in which the official IAU constellation boundaries are defined.
  * RA/Dec are converted to decimal values manually (not via a sexagesimal string
    parser), so slightly non-standard values such as arcseconds = 60.0 still work.
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import numpy as np
    import astropy.units as u
    from astropy.coordinates import SkyCoord
except ImportError:
    sys.exit("Missing dependencies. Install them with:  pip install astropy numpy")

# HHMMSS.ss(+|-)DDMMSS.s
COORD_RE = re.compile(
    r"^(?P<h>\d{2})(?P<m>\d{2})(?P<s>\d{2}(?:\.\d+)?)"
    r"(?P<sign>[+-])"
    r"(?P<d>\d{2})(?P<am>\d{2})(?P<as_>\d{2}(?:\.\d+)?)$"
)


def parse_coordinate(text):
    """Return (ra_hours, dec_degrees) as floats, or None if the string is invalid."""
    match = COORD_RE.match(text)
    if not match:
        return None
    g = match.groupdict()
    ra_h = int(g["h"]) + int(g["m"]) / 60.0 + float(g["s"]) / 3600.0
    dec = int(g["d"]) + int(g["am"]) / 60.0 + float(g["as_"]) / 3600.0
    if g["sign"] == "-":
        dec = -dec
    if not (0.0 <= ra_h < 24.0 and -90.0 <= dec <= 90.0):
        return None
    return ra_h, dec


def prompt_for_input_path():
    """Ask the user for the source file path until a valid file is given."""
    while True:
        try:
            raw = input("Enter the path of the source file: ").strip()
        except (EOFError, KeyboardInterrupt):
            sys.exit("\nNo input file provided. Exiting.")
        # Strip surrounding quotes (e.g. from drag-and-drop or 'Copy as path') and expand ~
        raw = raw.strip('"').strip("'")
        path = Path(raw).expanduser()
        if path.is_file():
            return path
        print(f"  File not found: {path}\n  Please try again (Ctrl+C to quit).")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("input", type=Path, nargs="?",
                        help="input text file (optional; you will be prompted if omitted)")
    parser.add_argument("output", type=Path, nargs="?", help="output file (default: <input>_out.txt)")
    parser.add_argument("--abbr", action="store_true",
                        help="write 3-letter IAU abbreviations instead of full names")
    args = parser.parse_args()
    if args.input is None:
        args.input = prompt_for_input_path()
    out_path = args.output or args.input.with_name(args.input.stem + "_out.txt")

    # --- Read input ---------------------------------------------------------
    try:
        with open(args.input, "r", encoding="utf-8") as fh:
            lines = [ln.strip() for ln in fh]
    except FileNotFoundError:
        sys.exit(f"Error: input file not found: {args.input}")
    except (PermissionError, UnicodeDecodeError, OSError) as exc:
        sys.exit(f"Error reading {args.input}: {exc}")

    # --- Parse (keep original line order; remember invalid lines) -----------
    entries = []          # (original_string, parsed_or_None) for non-empty lines
    for lineno, line in enumerate(lines, start=1):
        if not line:
            continue
        parsed = parse_coordinate(line)
        if parsed is None:
            print(f"Warning: line {lineno} could not be parsed: {line!r}", file=sys.stderr)
        entries.append((line, parsed))

    valid = [(i, p) for i, (_, p) in enumerate(entries) if p is not None]
    names = {}
    if valid:
        # --- Vectorised constellation lookup (fast for large files) ---------
        ra = np.array([p[0] for _, p in valid]) * u.hourangle
        dec = np.array([p[1] for _, p in valid]) * u.deg
        coords = SkyCoord(ra=ra, dec=dec, frame="icrs")
        result = coords.get_constellation(short_name=args.abbr)
        names = {idx: str(name) for (idx, _), name in zip(valid, np.atleast_1d(result))}

    # --- Write output -------------------------------------------------------
    try:
        with open(out_path, "w", encoding="utf-8") as fh:
            for i, (text, parsed) in enumerate(entries):
                label = names[i] if parsed is not None else "INVALID"
                fh.write(f"{text} ; {label}\n")
    except OSError as exc:
        sys.exit(f"Error writing {out_path}: {exc}")

    print(f"Processed {len(entries)} lines ({len(valid)} valid). Output: {out_path}")


if __name__ == "__main__":
    main()