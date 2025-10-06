import csv
import fnmatch
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


def discover_files(root_directory: str, patterns: Sequence[str]) -> List[str]:
    """
    Recursively discover files under the provided directory matching any of the glob-like patterns.

    Parameters
    ----------
    root_directory : str
        Root directory to scan.
    patterns : Sequence[str]
        Filename patterns (e.g., "*.fits", "*.pbcor.fits").

    Returns
    -------
    List[str]
        Absolute file paths.
    """
    normalized_patterns = [p.strip() for p in patterns if p and p.strip()]
    discovered: List[str] = []
    for current_root, _dirs, files in os.walk(root_directory):
        for file_name in files:
            for pattern in normalized_patterns:
                if fnmatch.fnmatch(file_name, pattern):
                    discovered.append(os.path.abspath(os.path.join(current_root, file_name)))
                    break
    # De-duplicate while preserving order
    seen: set = set()
    unique_files: List[str] = []
    for f in discovered:
        if f not in seen:
            unique_files.append(f)
            seen.add(f)
    return unique_files


# No FITS header access needed in this scanner per requirements.


def write_csv(
    rows: Iterable[Dict[str, Optional[object]]],
    output_csv: str,
    field_order: Sequence[str]
) -> None:
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(field_order))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in field_order})


def main() -> int:
    argv = sys.argv[1:]
    if len(argv) < 2:
        print("Usage: python scan_fits.py <directory> <output_csv> [patterns_csv]")
        print("Example: python scan_fits.py /srv/data out.csv '*cont*.pbcor.fits'\n")
        return 1

    directory = argv[0]
    output_csv = argv[1]
    patterns_csv = argv[2] if len(argv) >= 3 else "*.fits,*.FITS,*.pbcor.fits"

    patterns = [p.strip() for p in patterns_csv.split(",") if p.strip()]
    pattern_lower = patterns_csv.lower()
    if "cont" in pattern_lower:
        image_type = "cont"
    elif ("cube" in pattern_lower) or ("spw" in pattern_lower):
        image_type = "cube"
    else:
        image_type = "unknown"

    files = discover_files(directory, patterns)
    print(f"Discovered {len(files)} files matching patterns {patterns}")
    if not files:
        write_csv([], output_csv, [
            "path", "filename", "image_type", "file_size_bytes",
        ])
        print(f"No files found after filtering. Wrote empty CSV to {output_csv}")
        return 0

    rows: List[Dict[str, Optional[object]]] = []
    for idx, file_path in enumerate(files, start=1):
        row: Dict[str, Optional[object]] = {}
        row["path"] = file_path
        row["filename"] = os.path.basename(file_path)
        try:
            row["file_size_bytes"] = os.path.getsize(file_path)
        except Exception:
            row["file_size_bytes"] = None
        row["image_type"] = image_type
        rows.append(row)
        if idx % 200 == 0:
            print(f"Processed {idx} / {len(files)} files...")

    field_order = ["path", "filename", "image_type", "file_size_bytes"]
    write_csv(rows, output_csv, field_order)
    print(f"Wrote {len(rows)} rows to {output_csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())


