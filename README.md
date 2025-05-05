# CGAnalysis

A Python package for analyzing coarse-grained molecular dynamics simulations, especially designed for CALVADOS coarse-grained models and similar systems.
Supports cluster analysis, contact maps, inter-chain contact networks, MSD calculations, structure factor analysis and slab density map picture.

Note: Due to the storage limitation, the example case can be downloaded in [this link](http://ug.link/DXP4800PLUS-7C5/filemgr/share-download/?id=04f0522372ce443f93a8aa6d1ecd0d1b)

![image](https://github.com/user-attachments/assets/97f6e761-5773-4e26-bebb-0c149b32ed4e)

## Installation

Clone this repository and install in editable mode:

````bash
git clone https://github.com/aib001/CGAnalysis.git #you can use ssh for a quicker fetch
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

    z, kde = density.analyze_slab_density(slab_top, slab_traj, bandwidth=0.1, last_n_frames=200, n_workers=10)


```

## Demo
![contact_network_00100](https://github.com/user-attachments/assets/9dd4e0f6-8f2d-4e6e-84b3-6517ea2f62d9)
![contact_network_00200](https://github.com/user-attachments/assets/c45a423c-1787-4a88-8cb1-8b8565231742)
![structure_factor_q](https://github.com/user-attachments/assets/646407be-1d7f-45a3-a9f1-031e2e391573)
![violin_00300](https://github.com/user-attachments/assets/a260b71d-1066-4c0f-bea9-7890315a8ff2)


