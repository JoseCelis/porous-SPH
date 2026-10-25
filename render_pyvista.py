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
    .venv/bin/python render_pyvista.py part_1 --no-fluid              # rocks only, with the box outline
    .venv/bin/python render_pyvista.py part_1 --no-fluid --out rocks.png   # one still image (last frame)
"""
import argparse
import glob
import itertools
import os
import re
import sys

import numpy as np
import pyvista as pv

DIM = 0.63    # box half-width in x/y, matches plt.gnu
ZTOP = 2 * DIM  # z ceiling: sph.cpp puts it at zheight = 2*dim


def find_rocks_file(part_dir, explicit):
    """rocks.txt is written by sph.cpp next to the frames. Look in part_dir
    first (move rocks.txt along with the .dat files when archiving a run),
    then in the current directory."""
    if explicit:
        return explicit
    here = os.path.join(part_dir, "rocks.txt")
    if os.path.exists(here):
        return here
    if os.path.exists("rocks.txt") and os.path.abspath(part_dir) != os.path.abspath("."):
        print(f"warning: no rocks.txt in {part_dir}, using ./rocks.txt -- it may belong to a different run",
              file=sys.stderr)
        return "rocks.txt"
    return None


def load_rocks(path):
    """Read rocks.txt as written by setup_rocks() in sph.cpp. Returns
    (rocks, ramp_steps, steps_per_frame), rocks being a list of
    (center, faces) with faces rows "nx ny nz d" relative to the centre."""
    rocks, ramp, per_frame = [], 0, 1000
    with open(path) as f:
        for line in f:
            w = line.split()
            if not w or w[0].startswith("#"):
                continue
            if w[0] == "ramp_steps":
                ramp = int(w[1])
            elif w[0] == "steps_per_frame":
                per_frame = int(w[1])
            elif w[0] == "rock":
                rocks.append((np.array([float(v) for v in w[1:4]]), []))
            else:
                rocks[-1][1].append([float(v) for v in w[:4]])
    return [(c, np.array(fc)) for c, fc in rocks], ramp, per_frame


def rock_mesh(center, faces):
    """Surface of the convex polyhedron {p : nhat_k.(p-center) <= d_k}, built
    the same way sph.cpp defines it: every corner is where three face planes
    meet and lies inside all the others, and the surface is their convex hull."""
    n = faces[:, :3] / np.linalg.norm(faces[:, :3], axis=1)[:, None]
    d = faces[:, 3]
    verts = []
    for i, j, k in itertools.combinations(range(len(d)), 3):
        A = n[[i, j, k]]
        if abs(np.linalg.det(A)) < 1e-9:
            continue
        p = np.linalg.solve(A, d[[i, j, k]])
        if np.all(n @ p - d <= 1e-9):
            verts.append(p)
    verts = np.unique(np.round(verts, 12), axis=0)
    hull = pv.PolyData(verts).delaunay_3d().extract_surface(algorithm="dataset_surface")
    return hull.translate(center, inplace=False)


def rock_scale(path, ramp_steps, steps_per_frame):
    """sph.cpp grows the rocks from a point over rock_ramp_steps; frame NNN.dat
    is written at step (NNN-1)*steps_per_frame."""
    if ramp_steps <= 0 or path is None:
        return 1.0
    m = re.match(r"(\d+)\.dat$", os.path.basename(path))
    if not m:
        return 1.0
    ts = (int(m.group(1)) - 1) * steps_per_frame
    return min(1.0, ts / ramp_steps)


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
    ap.add_argument("--color", choices=["speed", "phase"], default="phase",
                     help="colour points by speed magnitude or phase (default: speed)")
    ap.add_argument("--out", default=None,
                     help="output path (default: <part_dir>/animation.mp4). A .png path saves one still "
                          "image of the last frame instead of a movie")
    ap.add_argument("--framerate", type=int, default=24)
    ap.add_argument("--azimuth", type=float, default=60.0, help="starting camera azimuth (deg)")
    ap.add_argument("--elevation", type=float, default=-20.0, help="starting camera elevation (deg)")
    ap.add_argument("--rotate", action="store_true", help="slowly rotate the camera through the animation")
    ap.add_argument("--rotate-degrees", type=float, default=360.0, help="total azimuth sweep if --rotate is set")
    ap.add_argument("--point-size", type=float, default=10.0)
    ap.add_argument("--opacity", type=float, default=1.0, help="point opacity, 0-1 (default: 1.0, opaque)")
    ap.add_argument("--limit", type=int, default=None, help="only render the first N frames (for quick tests)")
    ap.add_argument("--rocks", default=None,
                     help="rocks.txt written by sph.cpp (default: <part_dir>/rocks.txt, then ./rocks.txt)")
    ap.add_argument("--no-rock", action="store_true", help="do not draw the rocks")
    ap.add_argument("--rock-color", default="saddlebrown")
    ap.add_argument("--rock-opacity", type=float, default=1.0, help="rock opacity, 0-1 (default: 1.0, opaque)")
    ap.add_argument("--no-fluid", action="store_true",
                     help="draw only the rocks, no particles (the frames still set the rock growth and the "
                          "rotation; with no frames at all, the rocks are drawn once at full size)")
    box = ap.add_mutually_exclusive_group()
    box.add_argument("--box", dest="box", action="store_true", default=None,
                     help="draw the box outline (default: on with --no-fluid, off otherwise)")
    box.add_argument("--no-box", dest="box", action="store_false", help="do not draw the box outline")
    ap.add_argument("--window-size", type=int, nargs=2, default=[640, 480],
                     help="both dims should be multiples of 16 to avoid ffmpeg resizing (default: 896x704)")
    args = ap.parse_args()

    if args.no_fluid and args.no_rock:
        print("--no-fluid with --no-rock leaves nothing to draw", file=sys.stderr)
        sys.exit(1)
    show_box = args.no_fluid if args.box is None else args.box

    files = sorted(glob.glob(os.path.join(args.part_dir, "*.dat")))
    if not files and not args.no_fluid:
        print(f"no .dat files found in {args.part_dir}", file=sys.stderr)
        sys.exit(1)
    if args.limit:
        files = files[: args.limit]

    out_path = args.out or os.path.join(args.part_dir, "animation.mp4")
    still = out_path.lower().endswith(".png")
    if still and files:
        files = files[-1:]  # a still image shows the last frame only

    pl = pv.Plotter(off_screen=True, window_size=args.window_size)
    pl.set_background("white")

    cloud = None
    if not args.no_fluid:
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
            cmap = ["blue", "black"]
            label = "phase"

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

    if show_box:
        # the simulation box, so the rocks have a frame of reference without the fluid
        pl.add_mesh(pv.Box(bounds=(-DIM, DIM, -DIM, DIM, 0.0, ZTOP)).outline(), color="gray", line_width=1)

    rocks_path = None if args.no_rock else find_rocks_file(args.part_dir, args.rocks)
    rock = None
    if rocks_path:
        rock_list, rock_ramp, rock_per_frame = load_rocks(rocks_path)
        if rock_list:
            rock = True
    if args.no_fluid and rock is None:
        print(f"--no-fluid: no rocks to draw (no rocks.txt in {args.part_dir}, or it lists no rocks; "
              "use --rocks PATH)", file=sys.stderr)
        sys.exit(1)
    if rock is not None:
        # all rocks in one mesh (and one edge mesh), each point remembering its
        # own rock's centre so every rock can grow about that centre
        solids, edges = [], []
        for c, fc in rock_list:
            solid = rock_mesh(c, fc)
            # the hull is triangulated; outline only the real polyhedron edges
            # (the creases), not the diagonals that split each flat face
            edge = solid.extract_feature_edges(feature_angle=1.0, boundary_edges=False,
                                               non_manifold_edges=False, manifold_edges=False)
            solid.point_data["center"] = np.tile(c, (solid.n_points, 1))
            edge.point_data["center"] = np.tile(c, (edge.n_points, 1))
            solids.append(solid)
            edges.append(edge)
        rock_full = pv.merge(solids)
        rock_edges_full = pv.merge(edges)
        rock_actor = pl.add_mesh(rock_full.copy(), color=args.rock_color, opacity=args.rock_opacity,
                                 smooth_shading=False, scalars=None)
        edge_actor = pl.add_mesh(rock_edges_full.copy(), color="black", line_width=2, scalars=None)
        rock_shown = rock_actor.mapper.dataset
        edges_shown = edge_actor.mapper.dataset

        def set_rock(path):
            # the solver scales each rock's planes about its centre, i.e. a
            # uniform scaling of each solid about that point
            s = rock_scale(path, rock_ramp, rock_per_frame)
            for full, shown in ((rock_full, rock_shown), (rock_edges_full, edges_shown)):
                c = full.point_data["center"]
                shown.points = c + s * (full.points - c)
            rock_actor.SetVisibility(s > 0.0)
            edge_actor.SetVisibility(s > 0.0)

        set_rock(files[0] if files else None)

    if args.opacity < 1.0 or args.rock_opacity < 1.0:
        # Plain alpha blending draws translucent spheres in an arbitrary
        # (non-depth-sorted) order, which reintroduces the wrong-particle-
        # on-top artifact this script exists to avoid. Depth peeling fixes
        # that at the cost of extra render passes.
        pl.enable_depth_peeling(number_of_peels=8)

    pl.camera_position = "iso"
    pl.camera.azimuth = args.azimuth
    pl.camera.elevation = args.elevation

    if not files:
        files = [None]  # --no-fluid without frames: one picture of the full-size rocks

    if not still:
        pl.open_movie(out_path, framerate=args.framerate)

    n = len(files)
    step = args.rotate_degrees / n if args.rotate else 0.0
    for i, f in enumerate(files):
        if cloud is not None:
            pos, speed, phase = load_frame(f)
            cloud.points = pos
        if rock is not None:
            set_rock(f)
        if args.rotate:
            pl.camera.azimuth = args.azimuth + step * i
        if still:
            pl.screenshot(out_path)
        else:
            pl.write_frame()
        if (i + 1) % 100 == 0 or i == n - 1:
            print(f"  {i + 1}/{n} frames", file=sys.stderr)

    pl.close()
    print(f"wrote {out_path} ({'still image' if still else f'{n} frames'})")


if __name__ == "__main__":
    main()
