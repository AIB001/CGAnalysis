import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.signal import savgol_filter
import os
from multiprocessing import Pool

def _process_frame(args):
    frame_index, top, traj, chain_indices, default_box = args

    u = mda.Universe(top, traj)
    u.trajectory[frame_index]

    box = u.trajectory.ts.dimensions[:3]
    if np.any(box < 1.0):
        box = default_box.copy()

    if np.all(box > 1.0):
        u.atoms.wrap(compound='residues')

    positions = u.atoms.positions
    centroids = []
    for indices in chain_indices:
        centroid = positions[indices].mean(axis=0)
        centroids.append(centroid)

    time = u.trajectory.ts.time / 1000.0  # ns
    return frame_index, time, np.array(centroids)

def compute_msd_and_diffusion(top, traj, n_chains=500, chain_length=163, n_frames=500, n_workers=4):
    """
    Compute the mean squared displacement (MSD) and estimate the diffusion coefficient.
    """
    u = mda.Universe(top, traj)
    chain_indices = [u.residues[i * chain_length:(i + 1) * chain_length].atoms.indices for i in range(n_chains)]
    default_box = np.array([2000.0, 2000.0, 2000.0])

    frame_indices = []
    for ts in u.trajectory[:n_frames]:
        frame_indices.append(ts.frame)

    args_list = [(frame, top, traj, chain_indices, default_box) for frame in frame_indices]

    results = []
    with Pool(processes=n_workers) as pool:
        for result in tqdm(pool.imap(_process_frame, args_list), total=len(args_list), desc="Computing MSD (parallel)"):
            results.append(result)

    results.sort()
    times = [r[1] for r in results]
    centroids = [r[2] for r in results]

    times = np.array(times)
    centroids = np.array(centroids)  # shape: (frames, chains, 3)

    disp = centroids - centroids[0]
    msd = np.mean(np.sum(disp ** 2, axis=2), axis=1)

    if len(msd) > 11:
        msd_smooth = savgol_filter(msd, window_length=11, polyorder=3)
    else:
        msd_smooth = msd

    # Diffusion coefficient from linear fit to last 50% of data
    fit_start = int(0.5 * len(times))
    slope, intercept = np.polyfit(times[fit_start:], msd[fit_start:], 1)
    D = slope / 6.0

    ## ---- Plot ----
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 20, 'font.family': 'Times New Roman'})

    plt.figure(figsize=(8, 6))
    plt.plot(times, msd, 'o-', color='gray', alpha=0.4, label="Original MSD")
    plt.plot(times, msd_smooth, '-', color='blue', linewidth=2.5, label="Smoothed MSD")
    plt.plot(times, slope * times + intercept, 'r--', linewidth=2, label=f"Fit D={D:.2e} nm²/ns")

    plt.xlabel("Time (ns)", fontname='Times New Roman')
    plt.ylabel("MSD (nm²)", fontname='Times New Roman')
    plt.title("Mean Squared Displacement", fontname='Times New Roman')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(frameon=True, framealpha=0.8)

    out_dir = os.path.dirname(traj)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "msd_diffusion.png"), dpi=300)
    plt.close()

    return times, msd, D
