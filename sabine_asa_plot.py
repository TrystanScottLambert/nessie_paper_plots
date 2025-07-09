"""
Making the plot for Sabines talk - Fixed version
"""

import datetime

from nessie import FlatCosmology, RedshiftCatalog
from nessie.helper_funcs import create_density_function

import numpy as np
import pandas as pd
import pyvista as pv

from PIL import Image
import glob

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


def spherical_to_cartesian(ra_deg, dec_deg, z):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    # Comoving distance (Mpc)
    dist = astro_cosmo.comoving_distance(z).to(u.Mpc).value
    x = dist * np.cos(dec) * np.cos(ra)
    y = dist * np.cos(dec) * np.sin(ra)
    z = dist * np.sin(dec)
    return np.column_stack((x, y, z))

cosmo = FlatCosmology(0.7, 0.3)
astro_cosmo = FlatLambdaCDM(H0= 70, Om0= 0.3)

shark_mock = '/Users/00115372/Desktop/mock_catalogs/offical_waves_mocks/v0.3.0/wide/waves_wide_gals.parquet'

shark_data_frame = pd.read_parquet(shark_mock)
shark_data_frame = shark_data_frame[shark_data_frame['zobs'] < 0.2]

redshift = np.array(shark_data_frame['zobs'])
ra = np.array(shark_data_frame['ra'])
dec = np.array(shark_data_frame['dec'])
ab_mag = np.array(shark_data_frame['total_ab_dust_Z_VISTA'])
vel_errs = np.ones(len(shark_data_frame)) * 50

frac_area = 1133.86/41253
func = create_density_function(redshift, len(redshift), frac_area, cosmo)

red_cat = RedshiftCatalog(shark_data_frame['ra'], shark_data_frame['dec'], shark_data_frame['zobs'], func, cosmo)
now = datetime.datetime.now()
red_cat.run_fof(0.05, 18)
later = datetime.datetime.now()
shark_data_frame['group_id'] = red_cat.group_ids
print(f"time taken (n = {len(shark_data_frame)}): {later - now}")

unique_group_ids, counts= np.unique(red_cat.group_ids[red_cat.group_ids !=-1], return_counts = True)
cut = np.where(counts>2)
print("THE NUMBER OF GROUPS IS: ", len(cut[0]))
unique_group_ids = unique_group_ids[cut]

group_galaxies = shark_data_frame[shark_data_frame['group_id'].isin(unique_group_ids)]
group_catalog = (
    group_galaxies.groupby('group_id')
    .agg(
        ra_mean=('ra', 'median'),
        dec_mean=('dec', 'median'),
        zobs_mean=('zobs', 'median'),
        n_galaxies=('group_id', 'count')
    )
    .reset_index()
)

# Galaxy coordinates
gal_coords = spherical_to_cartesian(
    shark_data_frame['ra'],
    shark_data_frame['dec'],
    shark_data_frame['zobs']
)

# Group coordinates
group_coords = spherical_to_cartesian(
    group_catalog['ra_mean'],
    group_catalog['dec_mean'],
    group_catalog['zobs_mean']
)

print(f"Galaxy coordinates range: x={gal_coords[:, 0].min():.1f} to {gal_coords[:, 0].max():.1f}")
print(f"Group coordinates range: x={group_coords[:, 0].min():.1f} to {group_coords[:, 0].max():.1f}")

n_galaxies = group_catalog['n_galaxies'].values
sizes = 1.0 + 10.0 * (n_galaxies - n_galaxies.min()) / (n_galaxies.max() - n_galaxies.min())

# ---- PyVista Plotter ----

plotter = pv.Plotter(off_screen=True, window_size=(1080, 720))  # Square format, higher res
plotter.set_background("white")

# Add galaxies - back to simple rendering
gal_points = pv.PolyData(gal_coords)
plotter.add_points(gal_points, color="black", point_size=0.5, opacity=0.4)

# Add group spheres with scaled sizes - back to simple rendering
group_poly = pv.PolyData(group_coords)
group_poly['scale'] = sizes
sphere = pv.Sphere(radius=1.0)
glyphs = group_poly.glyph(scale='scale', geom=sphere)
plotter.add_mesh(glyphs, color='red', opacity=0.8)

# ---- Critical Fix: Set up proper camera position ----
# Calculate the center and bounds of your data
all_coords = np.vstack([gal_coords, group_coords])
center = np.mean(all_coords, axis=0)
bounds = np.array([
    [all_coords[:, 0].min(), all_coords[:, 0].max()],
    [all_coords[:, 1].min(), all_coords[:, 1].max()],
    [all_coords[:, 2].min(), all_coords[:, 2].max()]
])

# Set camera to look at the center of your data
plotter.camera.position = center + np.array([0, 0, np.ptp(all_coords[:, 2]) * 2])
plotter.camera.focal_point = center
plotter.camera.up = [0, 1, 0]

# Reset the camera to fit all data
plotter.reset_camera()

print(f"Camera position: {plotter.camera.position}")
print(f"Camera focal point: {plotter.camera.focal_point}")

# ---- Setup high-DPI frame saving ----
import os
frame_dir = "gif_frames"
if not os.path.exists(frame_dir):
    os.makedirs(frame_dir)

frame_count = 0

# ---- First rotation (360 degrees) ----
print("Starting first rotation...")
n_frames = 180  # Doubled for half speed (2 degrees per frame)
#for i in range(n_frames):
#    # Set absolute azimuth angle instead of relative
#    angle = i * (360.0 / n_frames)  # 2 degrees per frame
#    plotter.camera.azimuth = angle
#    plotter.render()  # Force render before writing frame
#    plotter.screenshot(f"{frame_dir}/frame_{frame_count:04d}.png", scale=2)  # 2x DPI
#    frame_count += 1
#    if i % 45 == 0:  # Print progress every ~90 degrees
#        print(f"  Frame {i+1}/{n_frames} completed (angle: {angle:.1f}°)")

# ---- Smooth zoom-in ----
print("Starting zoom-in...")
for i in range(75):  # Double the frames for twice as much zoom
    plotter.camera.zoom(1.015)  # Same zoom rate, but more frames = more total zoom
    plotter.render()  # Force render before writing frame
    plotter.screenshot(f"{frame_dir}/frame_{frame_count:04d}.png", scale=2)  # 2x DPI
    frame_count += 1
    if i % 10 == 0:
        print(f"  Zoom frame {i+1}/50 completed")

# ---- Second rotation at new zoom ----
print("Starting second rotation...")
for i in range(n_frames):
    # Set absolute azimuth angle instead of relative
    angle = i * (360.0 / n_frames)  # 2 degrees per frame
    plotter.camera.azimuth = angle
    plotter.render()  # Force render before writing frame
    plotter.screenshot(f"{frame_dir}/frame_{frame_count:04d}.png", scale=2)  # 2x DPI
    frame_count += 1
    if i % 45 == 0:
        print(f"  Frame {i+1}/{n_frames} completed (angle: {angle:.1f}°)")

# ---- Second zoom-in (double zoom) ----
#print("Starting second zoom-in...")
#for i in range(50):  # Another 50 frames for double zoom
#    plotter.camera.zoom(1.015)  # Same zoom rate
#    plotter.render()  # Force render before writing frame
#    plotter.screenshot(f"{frame_dir}/frame_{frame_count:04d}.png", scale=2)  # 2x DPI
#    frame_count += 1
#    if i % 10 == 0:
#        print(f"  Second zoom frame {i+1}/50 completed")

# ---- Third rotation at final zoom ----
#print("Starting third rotation...")
#for i in range(n_frames):
#    # Set absolute azimuth angle instead of relative
#    angle = i * (360.0 / n_frames)  # 2 degrees per frame
#    plotter.camera.azimuth = angle
#    plotter.render()  # Force render before writing frame
#    plotter.screenshot(f"{frame_dir}/frame_{frame_count:04d}.png", scale=2)  # 2x DPI
#    frame_count += 1
#    if i % 45 == 0:
#        print(f"  Frame {i+1}/{n_frames} completed (angle: {angle:.1f}°)")

# ---- Finalize ----
print("Creating GIF from high-DPI frames...")
plotter.close()

# Create GIF from the high-DPI PNG frames

# Get all frame files in order
frame_files = sorted(glob.glob(f"{frame_dir}/frame_*.png"))
print(f"Found {len(frame_files)} frames")

# Load frames and create GIF
frames = []
for i, frame_file in enumerate(frame_files):
    if i % 50 == 0:
        print(f"Processing frame {i+1}/{len(frame_files)}")
    img = Image.open(frame_file)
    frames.append(img)

# Save as GIF
print("Saving GIF...")
frames[0].save(
    "cosmic_groups_hd.gif",
    save_all=True,
    append_images=frames[1:],
    duration=100,  # 100ms per frame = 10 FPS
    loop=0
)

print("High-DPI GIF saved as cosmic_groups_hd.gif")
print(f"Total frames: {180 + 50 + 180 + 50 + 180} = 640 frames")
print(f"Frame files saved in '{frame_dir}' directory")