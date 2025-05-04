import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter
import os
from multiprocessing import Pool

def _process_frame(args):
    """ Helper function for parallel processing """
    frame_index, top, traj, contact_threshold, n_chains, chain_length = args

    u = mda.Universe(top, traj)
    u.trajectory[frame_index]

    chain_indices = [u.residues[i*chain_length:(i+1)*chain_length].atoms.indices for i in range(n_chains)]
    chain_ids = np.zeros(u.atoms.n_atoms, dtype=int)
    for chain_idx, indices in enumerate(chain_indices):
        chain_ids[indices] = chain_idx

    box = u.trajectory.ts.dimensions[:3]
    default_box = np.array([2000.0, 2000.0, 2000.0])
    if np.any(box < 1.0):
        box = default_box.copy()

    if np.all(box > 1.0):
        u.atoms.wrap(compound='residues')

    positions = u.atoms.positions

    if np.any(positions > box.max() + 1.0):
        raise ValueError(f"Frame {frame_index}: Atom positions exceed box after wrapping.")

    tree = cKDTree(positions, boxsize=box)
    pairs = tree.query_pairs(r=contact_threshold)

    count = 0
    for i, j in pairs:
        ci, cj = chain_ids[i], chain_ids[j]
        if ci != cj:
            count += 1

    return frame_index, count

def compute_interchain_contacts(top, traj, contact_threshold=8.0, n_chains=500, chain_length=163,
                                n_frames=500, process_interval=50, n_workers=4):
    """
    Compute the number of inter-chain contacts over time and plot the results (parallel version).
    """

    u = mda.Universe(top, traj)
    available_frames = [ts.frame for ts in u.trajectory[:n_frames] if ts.frame % process_interval == 0]

    # Prepare arguments for each frame
    frame_indices = []
    for ts in u.trajectory[:n_frames]:
        if ts.frame % process_interval == 0:
            frame_indices.append(ts.frame)

    args_list = [(frame, top, traj, contact_threshold, n_chains, chain_length) for frame in frame_indices]

    results = []
    with Pool(processes=n_workers) as pool:
        for result in tqdm(pool.imap(_process_frame, args_list), total=len(args_list), desc="Computing inter-chain contacts (parallel)"):
            results.append(result)

    results.sort()  # Sort by frame index
    frames, contact_counts = zip(*results)

    ## ------------------ Plotting ------------------

    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 20, 'font.family': 'Times New Roman'})

    time_ns = np.array(frames) * (u.trajectory.dt / 1000.0)  # convert to ns
    contacts = np.array(contact_counts)

    if len(contacts) > 11:
        contacts_smooth = savgol_filter(contacts, window_length=11, polyorder=3)
    else:
        contacts_smooth = contacts

    plt.figure(figsize=(10, 6))
    plt.plot(time_ns, contacts, 'o-', color='gray', alpha=0.4, label="Original")
    plt.plot(time_ns, contacts_smooth, '-', color='blue', linewidth=2.5, label="Smoothed")

    plt.xlabel("Time (ns)", fontname='Times New Roman')
    plt.ylabel("Inter-chain Contacts", fontname='Times New Roman')
    plt.title("Inter-chain Contact Evolution", fontname='Times New Roman')
    plt.legend()
    plt.tight_layout()

    out_dir = os.path.join(os.path.dirname(traj), "interchain_contact")
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, "interchain_contacts.png"), dpi=300)
    plt.close()

    return frames, contact_counts
