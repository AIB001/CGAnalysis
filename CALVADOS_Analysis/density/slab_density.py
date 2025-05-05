import MDAnalysis as mda
import numpy as np
from scipy.ndimage import uniform_filter1d
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm

def _process_frame_hist(args):
    """
    Process a single frame to compute the histogram of z-coordinates.
    """
    frame, top_path, dcd_path, bin_edges = args
    u = mda.Universe(top_path, dcd_path)
    u.trajectory[frame]
    ts = u.trajectory.ts

    z_coords = ts.positions[:, 2]
    hist, _ = np.histogram(z_coords, bins=bin_edges)

    return hist, z_coords, ts.dimensions[2]

def analyze_slab_density(top, traj, bin_width=5.0, last_n_frames=200, n_workers=4):
    """
    Analyze the slab density and compute the average density profile along z-axis.

    Returns:
    - z_centers: bin centers.
    - density_smooth: smoothed averaged density profile.
    """
    u = mda.Universe(top, traj)
    total_frames = len(u.trajectory)
    start_frame = max(total_frames - last_n_frames, 0)
    frames_to_process = range(start_frame, total_frames)

    u.trajectory[-1]
    z_min = u.atoms.positions[:, 2].min() - 5
    z_max = u.atoms.positions[:, 2].max() + 5
    bin_edges = np.arange(z_min, z_max + bin_width, bin_width)
    z_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    args_list = [(f, top, traj, bin_edges) for f in frames_to_process]

    histograms = []
    all_z_coords_list = []
    box_lengths = []

    with mp.Pool(processes=n_workers) as pool:
        for result in tqdm(pool.imap(_process_frame_hist, args_list), total=len(args_list), desc="Processing frames"):
            hist, z_coords, box_z = result
            histograms.append(hist)
            all_z_coords_list.append(z_coords)
            box_lengths.append(box_z)

    histograms = np.array(histograms)
    density_mean = np.mean(histograms, axis=0)
    density_mean = density_mean / (bin_width * len(all_z_coords_list[0]))

    density_smooth = uniform_filter1d(density_mean, size=2)

    # ---------- Auto-detect dense phase boundaries ----------
    first_deriv = np.gradient(density_smooth, z_centers)
    second_deriv = np.gradient(first_deriv, z_centers)
    zero_crossings = np.where(np.diff(np.sign(second_deriv)))[0]

    if len(zero_crossings) >= 2:
        first_deriv_abs = np.abs(first_deriv[zero_crossings])
        top_two_indices = np.argsort(first_deriv_abs)[-2:]
        top_two_indices_sorted = sorted(top_two_indices)
        z_start, z_end = z_centers[zero_crossings[top_two_indices_sorted]]
    else:
        z_start, z_end = z_centers[0], z_centers[-1]  # fallback to whole box

    avg_box_z = np.mean(box_lengths)
    length_ratio = (z_end - z_start) / avg_box_z if avg_box_z > 0 else 0.0

    all_z_coords = np.concatenate(all_z_coords_list)
    particle_ratio = np.sum((all_z_coords >= z_start) & (all_z_coords <= z_end)) / len(all_z_coords)

    print(f"Dense phase length ratio (average over last frames): {length_ratio:.4f}")
    print(f"Dense phase particle ratio (average over last frames): {particle_ratio:.4f}")

    # ---------- Plot ----------
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

    plt.figure(figsize=(14, 6))
    for hist in histograms:
        density = hist / (bin_width * len(all_z_coords_list[0]))
        plt.plot(z_centers, density, color='gray', alpha=0.3)

    plt.plot(z_centers, density_smooth, color='blue', linewidth=2.5, label='Average Density')

    plt.axvline(z_start, color='red', linestyle='--', label='Dense phase boundary')
    plt.axvline(z_end, color='red', linestyle='--')

    plt.xlabel('Z coordinate (Å)', fontname='Times New Roman')
    plt.ylabel('Density (arb. units)', fontname='Times New Roman')
    plt.title('Density profile along Z-axis', fontname='Times New Roman')
    plt.legend()
    plt.tight_layout()

    out_dir = os.path.dirname(traj)
    plt.savefig(os.path.join(out_dir, "z_density_profile.png"), dpi=300)
    plt.close()

    return z_centers, density_smooth
