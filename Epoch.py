from datetime import datetime, timezone


def iso_to_decimal_epoch(iso_ts: str) -> float:
    # Remove 'Z' if present
    iso_ts = iso_ts.rstrip("Z")

    # Truncate to microseconds and parsec
    dt = datetime.fromisoformat(iso_ts[:26])
    dt = dt.replace(tzinfo=timezone.utc)

    # Extract year and check for leap year
    year = dt.year
    is_leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
    days_in_year = 366 if is_leap else 365

    # Get day of year and fractional part of day
    doy = dt.timetuple().tm_yday
    seconds_in_day = dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6
    fractional_day = seconds_in_day / 86400

    # Compute decimal epoch
    decimal_epoch = year + (doy - 1 + fractional_day) / days_in_year
    return round(decimal_epoch, 6)


# Prompt user input
timestamp = input("timestamp: ")
epoch = iso_to_decimal_epoch(timestamp)
print("epoch:", epoch)
