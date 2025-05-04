import mdtraj as md
import numpy as np
import pandas as pd
import os
from tqdm import tqdm

def compute_contact_frequency(topol, traj, resname="LIG", distance_cutoff=0.4, n_frames=500):
    """
    Calculate contact frequency between ligand heavy atoms and protein residues.

    Parameters:
    - topol: Topology file path.
    - traj: Trajectory file path.
    - resname: Ligand residue name.
    - distance_cutoff: Contact distance cutoff (nm).
    - n_frames: Number of frames to analyze.

    Returns:
    - DataFrame with ligand atom indices and contact frequencies.
    """
    t = md.load(traj, top=topol)
    n_frames = min(n_frames, t.n_frames)

    ligand_atoms = t.topology.select(f"resname {resname} and element != H")
    protein_atoms = t.topology.select("protein and element != H")

    if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
        raise ValueError("Ligand or protein heavy atoms not found.")

    contact_counts = np.zeros(len(ligand_atoms), dtype=int)

    for frame in tqdm(range(n_frames), desc="Computing contacts"):
        distances = md.compute_distances(
            t.slice(frame),
            atom_pairs=[(i, j) for i in ligand_atoms for j in protein_atoms]
        )[0]

        for idx, dist in enumerate(distances):
            if dist < distance_cutoff:
                lig_atom = idx // len(protein_atoms)
                contact_counts[lig_atom] += 1

    frequency = contact_counts / n_frames
    result = pd.DataFrame({"Ligand_Atom": ligand_atoms, "Contact_Frequency": frequency})

    out_dir = os.path.dirname(os.path.abspath(traj))
    result.to_csv(os.path.join(out_dir, "ligand_contact_frequency.csv"), index=False)

    return result
