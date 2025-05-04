import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm

def compute_largest_size(top, traj, contact_threshold=6.0, n_chains=500, chain_length=163,
                         n_frames=500, process_interval=50):
    """
    Compute the size of the largest cluster at each frame.

    Parameters:
    - top: Topology file path.
    - traj: Trajectory file path.
    - contact_threshold: Contact distance threshold (Å).
    - n_chains: Number of chains.
    - chain_length: Number of residues per chain.
    - n_frames: Number of frames to analyze.
    - process_interval: Frame processing interval.
    """
    output_dir = os.path.join(os.path.dirname(traj), "largest_cluster_size")
    os.makedirs(output_dir, exist_ok=True)

    u = mda.Universe(top, traj)
    chain_indices = [u.residues[i * chain_length:(i + 1) * chain_length].atoms.indices for i in range(n_chains)]
    chain_ids = np.zeros(u.atoms.n_atoms, dtype=int)
    for chain_idx, indices in enumerate(chain_indices):
        chain_ids[indices] = chain_idx

    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 20, 'font.family': 'Times New Roman'})

    frames = []
    largest_sizes = []

    default_box = np.array([2000.0, 2000.0, 2000.0])  # 200 nm

    total_frames = len(u.trajectory[:n_frames])

    for ts in tqdm(u.trajectory[:n_frames], desc="Analyzing largest cluster size"):
        if ts.frame % process_interval != 0:
            continue

        box = ts.dimensions[:3]

        if np.any(box < 1.0):
            print(f"Warning: Frame {ts.frame} has invalid box dimensions {box}, using default box.")
            box = default_box.copy()

        if np.all(box > 1.0):
            u.atoms.wrap(compound='residues')

        positions = u.atoms.positions

        if np.any(positions > box.max() + 1.0):
            raise ValueError(f"Frame {ts.frame}: Atom positions exceed the box even after wrapping.")

        tree = cKDTree(positions, boxsize=box)
        pairs = tree.query_pairs(r=contact_threshold)

        adjacency = np.zeros((n_chains, n_chains), dtype=int)
        for i, j in pairs:
            ci, cj = chain_ids[i], chain_ids[j]
            if ci != cj:
                adjacency[ci, cj] = 1
                adjacency[cj, ci] = 1

        _, labels = connected_components(csr_matrix(adjacency), directed=False)
        cluster_sizes = np.bincount(labels)
        largest_size = cluster_sizes.max() if len(cluster_sizes) > 0 else 0

        frames.append(ts.frame)
        largest_sizes.append(largest_size)

    plt.figure()
    plt.plot(frames, largest_sizes, marker='s')
    plt.xlabel("Frame", fontname='Times New Roman')
    plt.ylabel("Largest Cluster Size", fontname='Times New Roman')
    plt.title("Largest Cluster Over Time", fontname='Times New Roman')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "largest_cluster_size.png"))
    plt.close()
