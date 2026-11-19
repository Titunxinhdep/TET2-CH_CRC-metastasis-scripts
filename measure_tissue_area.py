"""
Measure tissue area (the stained region) on stitched whole-section images
exported from a Keyence BZ-X microscope.

The micrometres-per-pixel calibration is read from the XML metadata embedded in
each .tif file, so images exported at different sizes are scaled correctly.

Expected folder layout:

    adoptive transfer/          <- point the script at this folder
        control 6-2/  1.tif, 1-t.tif, 2.tif, 2-t.tif, ...
        KO1-1/        1.tif, 1-t.tif, ...
        ...

Only files whose name is a plain number are measured (1.tif ... 10.tif);
files such as 1-t.tif and any other file are ignored.
Only folders whose name begins with "control" or "KO" are read; change
FOLDER_PREFIXES if your naming differs.
All sub-folders belonging to one animal are read (-1, -2, -3 ...; ten slides each).
Images excluded because of a technical fault are declared in EXCLUDE and are
not measured.

Install once:  pip3 install numpy pillow scikit-image scipy
Run:           python3 measure_tissue_area.py "/path/to/adoptive transfer"

Output, written to a "results" folder:
    area_per_image.csv       one row per image
    missing_calibration.csv  images with no calibration in the file
    total_by_folder.csv      total area and image count for each folder
    overlay/                 images with the detected boundary drawn in yellow,
                             for visual inspection
"""
import sys, re, struct, csv
from pathlib import Path
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu, gaussian
from skimage import measure

Image.MAX_IMAGE_PIXELS = None
MIN_PIECE_UM2 = 100000                 # discard fragments and debris below 0.1 mm2 (units: um2)
FOLDER_PREFIXES = ("control", "ko")    # only read folders starting with these (case-insensitive)
# Images excluded because of a technical fault: (folder name, file name). Add rows as needed.
EXCLUDE = {
    ("control 2-2", "9.tif"),          # control 2, slide 19
}
IMAGE_NAME = re.compile(r"^(\d+)\.tiff?$", re.IGNORECASE)   # 1.tif ... 10.tif; excludes 1-t.tif and macOS ._ files


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def read_keyence_um_per_px(path):
    """Return the Calibration value (nm per pixel of the saved image) from the BZ-X XML metadata."""
    raw = Path(path).read_bytes()
    i = raw.find(b"<?xml")
    if i < 0:
        return None
    j = raw.find(b"</Data>", i)
    xml = raw[i:j].decode("utf8", "ignore")
    m = re.search(r'<Calibration Type="System.Double">(-?\d+)</Calibration>', xml)
    if not m:
        return None
    nm = struct.unpack("<d", int(m.group(1)).to_bytes(8, "little", signed=True))[0]
    return nm / 1000.0


def tissue_mask(rgb, um_per_px):
    """Segment tissue on the saturation channel: Gaussian smoothing, Otsu threshold, size filter."""
    sat = np.asarray(Image.fromarray(rgb).convert("HSV"))[..., 1].astype(float)
    sat = gaussian(sat, 2, preserve_range=True)
    mask = sat > threshold_otsu(sat)
    min_px = max(1, int(MIN_PIECE_UM2 / um_per_px**2))
    lab = measure.label(mask)
    sizes = np.bincount(lab.ravel())
    keep = sizes >= min_px
    keep[0] = False
    return keep[lab]


def main(root):
    root = Path(root)
    out = root / "results"
    (out / "overlay").mkdir(parents=True, exist_ok=True)

    folders = sorted((d for d in root.iterdir()
                      if d.is_dir() and d.name.lower().startswith(FOLDER_PREFIXES)),
                     key=lambda d: natural_key(d.name))
    if not folders:
        sys.exit(f"No folder starting with {FOLDER_PREFIXES} was found in {root}")

    excluded = {(f.lower().replace(" ", ""), i.lower()) for f, i in EXCLUDE}
    rows, summary, no_calib = [], [], []
    for d in folders:
        images = sorted((p for p in d.iterdir() if IMAGE_NAME.match(p.name)),
                        key=lambda p: int(IMAGE_NAME.match(p.name).group(1)))
        total, total_filled, n_ok = 0.0, 0.0, 0
        for p in images:
            if (d.name.lower().replace(" ", ""), p.name.lower()) in excluded:
                print(f"  [EXCLUDED] {d.name}/{p.name}: listed in EXCLUDE")
                continue
            um = read_keyence_um_per_px(p)
            if um is None:
                print(f"  [SKIPPED] {d.name}/{p.name}: no calibration found in the file")
                no_calib.append({"folder": d.name, "image": p.name})
                continue
            rgb = np.asarray(Image.open(p).convert("RGB"))
            mask = tissue_mask(rgb, um)
            filled = ndi.binary_fill_holes(mask)
            px_mm2 = um**2 / 1e6
            area, area_f = mask.sum() * px_mm2, filled.sum() * px_mm2
            rows.append({"folder": d.name, "image": p.name, "um_per_px": round(um, 4),
                         "tissue_mm2": round(area, 3), "tissue_filled_mm2": round(area_f, 3),
                         "n_pieces": int(measure.label(mask).max())})
            total += area; total_filled += area_f; n_ok += 1

            edge = mask ^ ndi.binary_erosion(mask, iterations=3)
            ov = rgb.copy(); ov[edge] = [255, 255, 0]
            im = Image.fromarray(ov); im.thumbnail((1500, 1500))
            im.save(out / "overlay" / f"{d.name}_{p.stem}.jpg", quality=85)

        summary.append({"folder": d.name, "n_images": n_ok,
                        "total_tissue_mm2": round(total, 3),
                        "total_tissue_filled_mm2": round(total_filled, 3)})
        print(f"{d.name}: {n_ok} images, total {total:.2f} mm2")

    with open(out / "area_per_image.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["folder"])
        w.writeheader(); w.writerows(rows)
    with open(out / "total_by_folder.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader(); w.writerows(summary)
    with open(out / "missing_calibration.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["folder", "image"])
        w.writeheader(); w.writerows(no_calib)
    print(f"\nDone: {len(rows)} images in {len(folders)} folders -> {out}")
    if no_calib:
        print(f"{len(no_calib)} images had no calibration -> see missing_calibration.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
