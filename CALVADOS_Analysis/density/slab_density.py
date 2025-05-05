import MDAnalysis as mda
import numpy as np
from scipy.stats import gaussian_kde
import multiprocessing as mp
import matplotlib.pyplot as plt
import seaborn as sns
import os

def _process_frame(args):
    """
    Process a single frame to compute the dense phase proportion 
    and KDE curve along the z-axis.
    """
    frame, top_path, dcd_path, bandwidth = args
    u = mda.Universe(top_path, dcd_path)
    u.trajectory[frame]
    ts = u.trajectory.ts

    z_coords = ts.positions[:, 2]

    kde = gaussian_kde(z_coords, bw_method=bandwidth)
    z_min, z_max = z_coords.min(), z_coords.max()
    z_grid = np.linspace(z_min - 5, z_max + 5, 1000)
    density = kde.evaluate(z_grid)

    first_deriv = np.gradient(density, z_grid)
    second_deriv = np.gradient(first_deriv, z_grid)

    zero_crossings = np.where(np.diff(np.sign(second_deriv)))[0]

    if len(zero_crossings) >= 2:
        first_deriv_abs = np.abs(first_deriv[zero_crossings])
        top_two_indices = np.argsort(first_deriv_abs)[-2:]
        top_two_indices_sorted = sorted(top_two_indices)
        z_start, z_end = z_grid[zero_crossings[top_two_indices_sorted]]
    else:
        z_start, z_end = 0.0, 0.0

    total_z_length = ts.dimensions[2]
    phase_length = z_end - z_start
    length_ratio = phase_length / total_z_length if total_z_length != 0 else 0.0

    in_dense = (z_coords >= z_start) & (z_coords <= z_end)
    particle_ratio = np.sum(in_dense) / len(z_coords)

    return (length_ratio, particle_ratio, z_grid, density)

def analyze_slab_density(top, traj, bandwidth=0.1, last_n_frames=200, n_workers=4):
    """
    Analyze the slab density and compute the average dense phase proportions and KDE curves.
    
    Parameters:
    - top: Topology file path.
    - traj: Trajectory file path.
    - bandwidth: KDE bandwidth parameter.
    - last_n_frames: Number of frames from the end to analyze.
    - n_workers: Number of parallel processes.

    Returns:
    - final_z: z-grid array.
    - final_kde: averaged KDE curve over all processed frames.
    """
    u = mda.Universe(top, traj)
    total_frames = len(u.trajectory)
    start_frame = max(total_frames - last_n_frames, 0)
    frames_to_process = range(start_frame, total_frames)

    args_list = [(f, top, traj, bandwidth) for f in frames_to_process]

    with mp.Pool(processes=n_workers) as pool:
        results = pool.map(_process_frame, args_list)

    length_ratios, particle_ratios, z_arrays, kde_arrays = zip(*results)

    avg_length_ratio = np.mean(length_ratios)
    avg_particle_ratio = np.mean(particle_ratios)

    print(f"Dense phase length ratio (average over last frames): {avg_length_ratio:.4f}")
    print(f"Dense phase particle ratio (average over last frames): {avg_particle_ratio:.4f}")

    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

    plt.figure(figsize=(14, 6))
    for z, kde in zip(z_arrays, kde_arrays):
        plt.plot(z, kde, color='gray', alpha=0.3)

    z_common = z_arrays[0]
    kde_mean = np.mean(kde_arrays, axis=0)
    plt.plot(z_common, kde_mean, color='blue', linewidth=2.5, label='Average KDE')

    plt.xlabel('Z coordinate (Å)', fontname='Times New Roman')
    plt.ylabel('Density (arb. units)', fontname='Times New Roman')
    plt.title('KDE of particle positions along Z-axis', fontname='Times New Roman')
    plt.legend()
    plt.tight_layout()

    out_dir = os.path.dirname(traj)
    plt.savefig(os.path.join(out_dir, "z_density_kde.png"), dpi=300)
    plt.close()

    return z_common, kde_mean
