# CGAnalysis V1.0--By A.I.B. 

A Python package for analyzing coarse-grained molecular dynamics simulations, especially designed for CALVADOS coarse-grained models and similar systems.
Supports cluster analysis, contact maps, inter-chain contact networks, MSD calculations, structure factor analysis, slab density map picture as well as Flory-Huggins Theory Fitting.

Note: Due to the storage limitation, the example case can be downloaded in [this link](http://ug.link/DXP4800PLUS-7C5/filemgr/share-download/?id=04f0522372ce443f93a8aa6d1ecd0d1b)
![image](https://github.com/user-attachments/assets/999775c5-84fe-414c-af0f-ea1e3dbb6eec)


![image](https://github.com/user-attachments/assets/97f6e761-5773-4e26-bebb-0c149b32ed4e)

## Installation

Clone this repository and install in editable mode:

````bash
git clone https://github.com/aib001/CGAnalysis.git #you can use ssh for a quicker download
cd CGAnalysis
pip install -e .
````

Dependencies (listed in `requirements.txt`):

- MDAnalysis
- numpy
- scipy
- matplotlib
- seaborn
- tqdm
- networkx
- sklearn

To install:

````bash
pip install -r requirements.txt
````

## Quick Start

```python
if __name__ == '__main__': # To ensure multi_process begin successfully, the code should run under '__main__'
    from CALVADOS_Analysis.cluster import size_distribution, time_evolution, largest_size
    from CALVADOS_Analysis.contact import frequency, heatmap, interchain
    from CALVADOS_Analysis.diffusion import compute_msd_and_diffusion
    from CALVADOS_Analysis.structure_factor import compute_structure_factor

    ##############################
    # Define the path
    ##############################
    import os
    base_dir = os.path.dirname(os.path.abspath(__file__))
    top = os.path.join(base_dir, "FUS_1_163_500", "checkpoint.pdb")
    traj = os.path.join(base_dir, "FUS_1_163_500", "FUS_1_163_500_quench_quench.dcd")

    n_frames = 2000

    ##############################
    # cluster Analysis
    ##############################
    size_distribution.analyze_clusters(top, traj, n_frames=n_frames)
    time_evolution.analyze_large_clusters(top, traj, n_frames=n_frames)
    largest_size.compute_largest_size(top, traj, n_frames=n_frames)

    ##############################
    # Contact Analysis
    ##############################
    heatmap.generate_contact_heatmaps(top, traj, n_frames=n_frames)
    frames, contacts = interchain.compute_interchain_contacts(top, traj, n_frames=n_frames)

    ##############################
    # Structural Factor and MSD Analysis
    ##############################
    time_q, q_avg = compute_structure_factor(top, traj, n_frames=2000, n_windows=200, num_grid=32, n_workers=15)
    times, msd, D = compute_msd_and_diffusion(top, traj, n_frames=n_frames, n_workers=20)


    #############################
    # Slab Density Calculation
    #############################
    from CALVADOS_Analysis import density

    slab_top = os.path.join(base_dir, "slab_case", "checkpoint.pdb")
    slab_traj =os.path.join(base_dir, "slab_case", "FUS_1_163_500.dcd")

    z, kde = density.analyze_slab_density(slab_top, slab_traj, bandwidth=5.0, last_n_frames=200, n_workers=10)

    #############################
    # Flory-Huggins Theory Fitting
    #############################

    from CALVADOS_Analysis.density import flory_huggins

    T_all = [280, 285, 290, 295, 300, 305]
    rho_H = [494.29, 468.50, 443.09, 415.71, 373.94, 341.96]
    rho_L = [7.95, 9.39, 12.05, 11.89, 24.04, 26.68]

    flory_huggins.analyze_flory_huggins(temp=T_all, rho_H=rho_H, rho_L=rho_L, rho_protein=900)


```

## Demo
![z_density_profile](https://github.com/user-attachments/assets/552475fb-0599-4fc6-b2fd-25f642ba7070)
![contact_network_00100](https://github.com/user-attachments/assets/9dd4e0f6-8f2d-4e6e-84b3-6517ea2f62d9)
![contact_network_00200](https://github.com/user-attachments/assets/c45a423c-1787-4a88-8cb1-8b8565231742)
![structure_factor_q](https://github.com/user-attachments/assets/646407be-1d7f-45a3-a9f1-031e2e391573)
![violin_00300](https://github.com/user-attachments/assets/a260b71d-1066-4c0f-bea9-7890315a8ff2)
![largest_cluster_size](https://github.com/user-attachments/assets/1c9232a2-d89d-4055-8e99-f3895597df77)



