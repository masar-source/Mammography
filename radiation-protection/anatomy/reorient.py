"""Rebuild an included GATE label image from its original Radboud NRRD.

Usage: python reorient.py original/KBCT...nrrd output.nrrd
The native array convention is (i, j, k) with i fastest on disk.
"""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np


def convert(source, target):
    header = {}
    with Path(source).open("rb") as fp:
        if not fp.readline().startswith(b"NRRD"):
            raise ValueError("Expected a NRRD volume.")
        while line := fp.readline().strip():
            if line.startswith(b"#"):
                continue
            key, text = line.decode("ascii").split(": ", 1)
            header[key] = text
        if (header["type"], header["dimension"], header["encoding"]) != ("unsigned char", "3", "gzip"):
            raise ValueError("Expected an unsigned-byte 3D gzip NRRD from the CC archive.")
        vectors = [list(map(float, part.split(","))) for part in re.findall(r"\(([^()]*)\)", header["space directions"])]
        directions = np.array(vectors)
        spacing = np.diag(directions)
        if not (np.allclose(directions, np.diag(spacing)) and np.allclose(spacing, 0.273, atol=1e-7)):
            raise ValueError("The source orientation or voxel spacing differs from the verified release.")
        sizes = tuple(map(int, header["sizes"].split()))
        with gzip.GzipFile(fileobj=fp) as compressed:
            a = np.frombuffer(compressed.read(), dtype=np.uint8).reshape(sizes, order="F")
    if a.max() > 3:
        raise ValueError("Unrecognized tissue label.")
    counts = [int(np.count_nonzero(a == label)) for label in range(4)]
    mask = a != 0
    bounds = []
    for k in range(3):
        profile = mask.any(axis=tuple(i for i in range(3) if i != k))
        indices = np.flatnonzero(profile)
        bounds.append((int(indices[0]), int(indices[-1]) + 1))
    slab = a[tuple(slice(lo, hi) for lo, hi in bounds)]
    # Source +i = posterior to anterior and +j = upwards from the support.
    # Right-handed rotation into GATE: (+i, -k, +j). No interpolation.
    gate = slab[:, :, ::-1].transpose(0, 2, 1)
    ds = spacing[[0, 2, 1]]
    header_lines = [
        "NRRD0004", "type: unsigned char", "dimension: 3", "space dimension: 3",
        "sizes: " + " ".join(map(str, gate.shape)),
        "space directions: " + " ".join("(" + ",".join(map(str, v)) + ")" for v in np.diag(ds)),
        "kinds: domain domain domain", "encoding: gzip", "space origin: (0,0,0)", "", "",
    ]
    with Path(target).open("wb") as output:
        output.write("\n".join(header_lines).encode("ascii"))
        with gzip.GzipFile(filename="", fileobj=output, mode="wb", compresslevel=6, mtime=0) as zipped:
            for plane in gate.transpose(2, 1, 0):
                zipped.write(plane.tobytes())
    return dict(
        source_dimensions_ijk=list(sizes), source_space=header["space"],
        source_origin=header["space origin"], source_directions=header["space directions"],
        source_tissue_labels_voxels=counts, crop_bounds_ijk_half_open=bounds,
        gate_dimensions_xyz=list(gate.shape), gate_spacing_mm=ds.tolist(),
        tissue_extent_xyz_mm=(np.array(gate.shape)*ds).tolist(),
        total_tissue_cm3=(sum(counts[1:]) * float(np.prod(ds)) / 1000),
        glandular_volume_fraction=counts[2] / (counts[1] + counts[2]),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    args = parser.parse_args()
    if args.source.resolve() == args.target.resolve():
        parser.error("Choose a distinct output file.")
    convert(args.source, args.target)
