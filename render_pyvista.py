#!/usr/bin/env python3
"""Render a part_N/ directory of SPH .dat frames to an mp4.

Uses PyVista (VTK) so points get real z-buffered occlusion when the
camera is rotated -- unlike gnuplot's splot, which only respects plot
order.

Data format (written by kinetic_E() in sph.cpp), comma separated:
    x, y, z, vx, vy, vz, phase

Usage:
    .venv/bin/python render_pyvista.py part_1
    .venv/bin/python render_pyvista.py part_1 --color phase
    .venv/bin/python render_pyvista.py part_1 --azimuth 45 --elevation 20
    .venv/bin/python render_pyvista.py part_1 --rotate --limit 50 --out test.mp4
"""
import argparse
import glob
import os
import sys

import numpy as np
import pyvista as pv

DIM = 0.63    # box half-width in x/y, matches plt.gnu
ZTOP = 3.5    # z ceiling, matches plt.gnu


def load_frame(path):
    data = np.loadtxt(path, delimiter=",")
    pos = np.ascontiguousarray(data[:, 0:3])
    vel = data[:, 3:6]
    phase = data[:, 6]
    speed = np.linalg.norm(vel, axis=1)
    return pos, speed, phase


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("part_dir", help="directory containing NNN.dat frames, e.g. part_1")
    ap.add_argument("--color", choices=["speed", "phase"], default="speed",
                     help="colour points by speed magnitude or phase (default: speed)")
    ap.add_argument("--out", default=None, help="output mp4 path (default: <part_dir>/animation.mp4)")
    ap.add_argument("--framerate", type=int, default=24)
    ap.add_argument("--azimuth", type=float, default=45.0, help="starting camera azimuth (deg)")
    ap.add_argument("--elevation", type=float, default=20.0, help="starting camera elevation (deg)")
    ap.add_argument("--rotate", action="store_true", help="slowly rotate the camera through the animation")
    ap.add_argument("--rotate-degrees", type=float, default=360.0, help="total azimuth sweep if --rotate is set")
    ap.add_argument("--point-size", type=float, default=10.0)
    ap.add_argument("--opacity", type=float, default=1.0, help="point opacity, 0-1 (default: 1.0, opaque)")
    ap.add_argument("--limit", type=int, default=None, help="only render the first N frames (for quick tests)")
    ap.add_argument("--window-size", type=int, nargs=2, default=[640, 480],
                     help="both dims should be multiples of 16 to avoid ffmpeg resizing (default: 896x704)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.part_dir, "*.dat")))
    if not files:
        print(f"no .dat files found in {args.part_dir}", file=sys.stderr)
        sys.exit(1)
    if args.limit:
        files = files[: args.limit]

    out_path = args.out or os.path.join(args.part_dir, "animation.mp4")

    pos0, speed0, phase0 = load_frame(files[0])
    cloud = pv.PolyData(pos0)

    if args.color == "speed":
        cloud["scalar"] = speed0
        clim = [0, 6]
        cmap = "coolwarm"
        label = "|v|  (m/s)"
    else:
        cloud["scalar"] = phase0
        clim = [1, 2]
        cmap = ["blue", "yellow"]
        label = "phase"

    pl = pv.Plotter(off_screen=True, window_size=args.window_size)
    pl.set_background("white")
    pl.add_mesh(
        cloud,
        # scalars="scalar",
        cmap=cmap,
        # clim=clim,
        point_size=args.point_size,
        render_points_as_spheres=True,
        show_scalar_bar=False,
        opacity=args.opacity,
        # scalar_bar_args={"title": label, "color": "black"},
    )
    if args.opacity < 1.0:
        # Plain alpha blending draws translucent spheres in an arbitrary
        # (non-depth-sorted) order, which reintroduces the wrong-particle-
        # on-top artifact this script exists to avoid. Depth peeling fixes
        # that at the cost of extra render passes.
        pl.enable_depth_peeling(number_of_peels=8)

    pl.camera_position = "iso"
    pl.camera.azimuth = args.azimuth
    pl.camera.elevation = args.elevation

    pl.open_movie(out_path, framerate=args.framerate)

    n = len(files)
    step = args.rotate_degrees / n if args.rotate else 0.0
    for i, f in enumerate(files):
        pos, speed, phase = load_frame(f)
        cloud.points = pos
        if args.rotate:
            pl.camera.azimuth = args.azimuth + step * i
        pl.write_frame()
        if (i + 1) % 100 == 0 or i == n - 1:
            print(f"  {i + 1}/{n} frames", file=sys.stderr)

    pl.close()
    print(f"wrote {out_path} ({n} frames)")


if __name__ == "__main__":
    main()
