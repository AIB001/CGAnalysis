# import MDAnalysis as mda
# import numpy as np
# from scipy.spatial import cKDTree
# import matplotlib.pyplot as plt
# import seaborn as sns
# import os
# from tqdm import tqdm
# import networkx as nx
# from scipy.sparse import csr_matrix
# from scipy.sparse.csgraph import connected_components
# from matplotlib.colors import LinearSegmentedColormap
# from matplotlib.cm import ScalarMappable

# def generate_contact_heatmaps(top, traj, contact_threshold=6.0, n_chains=500, chain_length=163,
#                               n_frames=500, process_interval=50):
#     """
#     Generate contact heatmaps and contact frequency network graphs for chain-chain interactions.
#     """
#     output_dir = os.path.join(os.path.dirname(traj), "contact_heatmaps")
#     os.makedirs(output_dir, exist_ok=True)

#     network_dir = os.path.join(os.path.dirname(traj), "contact_networks")
#     os.makedirs(network_dir, exist_ok=True)

#     u = mda.Universe(top, traj)
#     chain_indices = [u.residues[i*chain_length:(i+1)*chain_length].atoms.indices for i in range(n_chains)]
#     chain_ids = np.zeros(u.atoms.n_atoms, dtype=int)
#     for chain_idx, indices in enumerate(chain_indices):
#         chain_ids[indices] = chain_idx

#     sns.set_style("whitegrid")
#     plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

#     default_box = np.array([2000.0, 2000.0, 2000.0])
#     chain_contact_frequency = np.zeros(n_chains, dtype=int)

#     cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])

#     for ts in tqdm(u.trajectory[:n_frames], desc="Generating contact heatmaps and networks"):
#         if ts.frame % process_interval != 0:
#             continue

#         box = ts.dimensions[:3]
#         if np.any(box < 1.0):
#             print(f"Warning: Frame {ts.frame} has invalid box {box}, using default box.")
#             box = default_box.copy()

#         if np.all(box > 1.0):
#             u.atoms.wrap(compound='residues')

#         positions = u.atoms.positions
#         if np.any(positions > box.max() + 1.0):
#             raise ValueError(f"Frame {ts.frame}: Atom positions exceed box after wrapping.")

#         tree = cKDTree(positions, boxsize=box)
#         pairs = tree.query_pairs(r=contact_threshold)

#         contact_matrix = np.zeros((n_chains, n_chains), dtype=int)
#         contacts_this_frame = {}

#         for i, j in pairs:
#             ci, cj = chain_ids[i], chain_ids[j]
#             if ci != cj:
#                 contact_matrix[ci, cj] += 1
#                 contact_matrix[cj, ci] += 1
#                 chain_contact_frequency[ci] += 1
#                 chain_contact_frequency[cj] += 1
#                 contacts_this_frame.setdefault(ci, set()).add(cj)
#                 contacts_this_frame.setdefault(cj, set()).add(ci)

#         ## ---- Heatmap ----
#         plt.figure(figsize=(10, 8))
#         sns.heatmap(contact_matrix, cmap="viridis", square=True,
#                     cbar_kws={'label': 'Contact count'})
#         plt.title(f"Contact Map - Frame {ts.frame}", fontname='Times New Roman')
#         plt.xlabel("Chain ID", fontname='Times New Roman')
#         plt.ylabel("Chain ID", fontname='Times New Roman')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir, f"contact_heatmap_{ts.frame:05d}.png"))
#         plt.close()

#         ## ---- Network ----
#         G = nx.Graph()
#         for ci in range(n_chains):
#             G.add_node(ci)

#         for ci in range(n_chains):
#             neighbors = contacts_this_frame.get(ci, set())
#             for cj in neighbors:
#                 if cj > ci:  # avoid double counting
#                     G.add_edge(ci, cj)

#         # --- Node positions based on real space centroids ---
#         centroids = []
#         for ci in range(n_chains):
#             atom_indices = chain_indices[ci]
#             centroids.append(u.atoms[atom_indices].center_of_mass())

#         centroids = np.array(centroids)
#         raw_pos = {ci: (centroids[ci, 0], centroids[ci, 1]) for ci in range(n_chains)}

#         # Normalize positions to 0~10 box for display
#         x_vals = centroids[:, 0]
#         y_vals = centroids[:, 1]
#         x_vals = (x_vals - x_vals.min()) / (x_vals.max() - x_vals.min())
#         y_vals = (y_vals - y_vals.min()) / (y_vals.max() - y_vals.min())
#         pos = {ci: (x_vals[ci]*10, y_vals[ci]*10) for ci in range(n_chains)}

#         ## ---- Fix long-distance edge by pulling connected nodes closer ----
#         max_allowed_distance = min(box[0], box[1]) / 2

#         for u_, v_ in G.edges():
#             # Calculate real wrapped distance
#             dvec = centroids[u_] - centroids[v_]
#             dvec -= box * np.round(dvec / box)  # minimum image convention
#             dist = np.linalg.norm(dvec)
#             if dist > max_allowed_distance:
#                 # Pull nodes closer in the plot
#                 pos[u_] = pos[v_]  # Force them overlap (could also average)

#         ## ---- Edge colors and weights ----
#         edge_colors = []
#         edge_weights = []

#         contact_values = [contact_matrix[u_, v_] for u_, v_ in G.edges()]
#         max_contact = max(contact_values) if contact_values else 1

#         for u_, v_ in G.edges():
#             strength = contact_matrix[u_, v_] / max_contact
#             edge_colors.append(cmap(strength))
#             edge_weights.append(0.5 + 2 * strength)

#         ## ---- Draw network ----
#         fig, ax = plt.subplots(figsize=(16, 16))
#         nodes = nx.draw_networkx_nodes(G, pos,
#                                        node_color="skyblue",
#                                        node_size=120,
#                                        ax=ax)
#         edges = nx.draw_networkx_edges(G, pos,
#                                        edge_color=edge_colors,
#                                        width=edge_weights,
#                                        alpha=0.8,
#                                        ax=ax)
#         labels = nx.draw_networkx_labels(G, pos,
#                                          font_size=8,
#                                          font_family='Times New Roman',
#                                          ax=ax)

#         ax.set_title(f"Contact Network - Frame {ts.frame}", fontname='Times New Roman')
#         ax.axis('off')

#         # Add colorbar
#         sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_contact))
#         sm.set_array([])
#         cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
#         cbar.set_label("Contact Strength", fontname='Times New Roman')

#         plt.tight_layout()
#         plt.savefig(os.path.join(network_dir, f"contact_network_{ts.frame:05d}.png"), dpi=300)
#         plt.close()


##########################
# Version V1.1
##########################

import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from multiprocessing import Pool

def _process_frame(args):
    frame_index, top, traj, contact_threshold, n_chains, chain_indices, default_box = args

    u = mda.Universe(top, traj)
    u.trajectory[frame_index]

    if np.all(u.trajectory.ts.dimensions[:3] > 1.0):
        u.atoms.wrap(compound='residues')

    box = u.trajectory.ts.dimensions[:3]
    if np.any(box < 1.0):
        box = default_box.copy()

    positions = u.atoms.positions

    # 质心 (NumPy 矢量化)
    centroids = np.zeros((n_chains, 3))
    for i, indices in enumerate(chain_indices):
        centroids[i] = positions[indices].mean(axis=0)

    # Contact matrix
    tree = cKDTree(positions, boxsize=box)
    pairs = tree.query_pairs(r=contact_threshold)

    contact_matrix = np.zeros((n_chains, n_chains), dtype=int)
    for i, j in pairs:
        ci, cj = 0, 0
        for idx, chain in enumerate(chain_indices):
            if i in chain:
                ci = idx
            if j in chain:
                cj = idx
        if ci != cj:
            contact_matrix[ci, cj] += 1
            contact_matrix[cj, ci] += 1

    return frame_index, contact_matrix, centroids, box

def generate_contact_heatmaps(top, traj, contact_threshold=6.0, n_chains=500, chain_length=163,
                              n_frames=500, process_interval=50, n_workers=10):
    """
    Parallelized: generate contact heatmaps and networks.
    """

    output_dir = os.path.join(os.path.dirname(traj), "contact_heatmaps")
    os.makedirs(output_dir, exist_ok=True)

    network_dir = os.path.join(os.path.dirname(traj), "contact_networks")
    os.makedirs(network_dir, exist_ok=True)

    u = mda.Universe(top, traj)
    chain_indices = [u.residues[i*chain_length:(i+1)*chain_length].atoms.indices for i in range(n_chains)]
    default_box = np.array([2000.0, 2000.0, 2000.0])

    frame_indices = []
    for ts in u.trajectory[:n_frames]:
        if ts.frame % process_interval == 0:
            frame_indices.append(ts.frame)

    args_list = [(frame, top, traj, contact_threshold, n_chains, chain_indices, default_box) for frame in frame_indices]

    # 并行计算每帧 contact_matrix 和 质心
    results = []
    with Pool(processes=n_workers) as pool:
        for result in tqdm(pool.imap(_process_frame, args_list), total=len(args_list),
                           desc="Calculating contacts and centroids (parallel)"):
            results.append(result)

    # colormap
    cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

    # 绘图（主线程）
    for frame_index, contact_matrix, centroids, box in tqdm(results, desc="Plotting heatmaps and networks"):

        ## ---- Heatmap ----
        plt.figure(figsize=(10, 8))
        sns.heatmap(contact_matrix, cmap="viridis", square=True,
                    cbar_kws={'label': 'Contact count'})
        plt.title(f"Contact Map - Frame {frame_index}", fontname='Times New Roman')
        plt.xlabel("Chain ID", fontname='Times New Roman')
        plt.ylabel("Chain ID", fontname='Times New Roman')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"contact_heatmap_{frame_index:05d}.png"))
        plt.close()

        ## ---- Network ----
        G = nx.Graph()
        for ci in range(n_chains):
            G.add_node(ci)

        for ci in range(n_chains):
            for cj in range(ci+1, n_chains):
                if contact_matrix[ci, cj] > 0:
                    G.add_edge(ci, cj)

        # Normalize positions
        x_vals = centroids[:, 0]
        y_vals = centroids[:, 1]
        x_vals = (x_vals - x_vals.min()) / (x_vals.max() - x_vals.min())
        y_vals = (y_vals - y_vals.min()) / (y_vals.max() - y_vals.min())
        pos = {ci: (x_vals[ci]*10, y_vals[ci]*10) for ci in range(n_chains)}

        ## ---- Fix distant pairs by PBC minimum image ----
        max_allowed_distance = min(box[0], box[1]) / 2
        for u_, v_ in G.edges():
            dvec = centroids[u_] - centroids[v_]
            dvec -= box * np.round(dvec / box)
            dist = np.linalg.norm(dvec)
            if dist > max_allowed_distance:
                pos[u_] = pos[v_]  # 拉近

        ## ---- Edge colors and weights ----
        edge_colors = []
        edge_weights = []

        contact_values = [contact_matrix[u_, v_] for u_, v_ in G.edges()]
        max_contact = max(contact_values) if contact_values else 1

        for u_, v_ in G.edges():
            strength = contact_matrix[u_, v_] / max_contact
            edge_colors.append(cmap(strength))
            edge_weights.append(0.5 + 2 * strength)

        ## ---- Draw network ----
        fig, ax = plt.subplots(figsize=(16, 16))
        nx.draw_networkx_nodes(G, pos, node_color="skyblue", node_size=120, ax=ax)
        nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_weights, alpha=0.8, ax=ax)
        nx.draw_networkx_labels(G, pos, font_size=8, font_family='Times New Roman', ax=ax)

        ax.set_title(f"Contact Network - Frame {frame_index}", fontname='Times New Roman')
        ax.axis('off')

        sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=max_contact))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
        cbar.set_label("Contact Strength", fontname='Times New Roman')

        plt.tight_layout()
        plt.savefig(os.path.join(network_dir, f"contact_network_{frame_index:05d}.png"), dpi=300)
        plt.close()
