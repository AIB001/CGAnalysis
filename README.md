# CGAnalysis

A Python package for analyzing coarse-grained molecular dynamics simulations, especially designed for CALVADOS coarse-grained models and similar systems.
Supports cluster analysis, contact maps, inter-chain contact networks, MSD calculations, and structure factor analysis.
Optimized for parallel computation and publication-quality visualization.

![Figure4_Phase_Diagram_01(1)](https://github.com/user-attachments/assets/aa764278-9d3a-4ed4-8a55-6bbb94c2a4e5)

## Package Structure

CGAnalysis/
 ├── CALVADOS_Analysis/
 │   ├── cluster/
 │   │   ├── size_distribution.py
 │   │   ├── time_evolution.py
 │   │   ├── largest_size.py
 │   ├── contact/
 │   │   ├── frequency.py
 │   │   ├── heatmap.py
 │   │   ├── interchain.py
 │   ├── diffusion/
 │   │   ├── msd_analysis.py
 │   ├── structure_factor/
 │   │   ├── structure_factor.py
 │   └── **init**.py
 ├── example/
 │   ├── example.py
 ├── setup.py
 ├── README.md
 ├── requirements.txt

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
from CALVADOS_Analysis.cluster import size_distribution, time_evolution, largest_size
from CALVADOS_Analysis.contact import frequency, heatmap, interchain
from CALVADOS_Analysis.diffusion import msd_analysis
from CALVADOS_Analysis.structure_factor import structure_factor

path = r'E:\My_Project\Phase_Seperation\FUS-PLD\FUS_1_163_2000_60cbbox\2000_60cbbob_h350_l300\FUS_1_163_2000'
top = f"{path}/checkpoint.pdb"
traj = f"{path}/FUS_1_163_2000_quench.dcd"

n_frames = 500

# Cluster analysis
size_distribution.analyze_clusters(top, traj, n_frames=n_frames)
time_evolution.analyze_large_clusters(top, traj, n_frames=n_frames)
largest_size.compute_largest_size(top, traj, n_frames=n_frames)

# Contact analysis
heatmap.generate_contact_heatmaps(top, traj, n_frames=n_frames, n_workers=4)
frames, contacts = interchain.compute_interchain_contacts(top, traj, n_frames=n_frames, n_workers=4)

# Diffusion & MSD
times, msd, D = msd_analysis.compute_msd_and_diffusion(top, traj, n_frames=n_frames, n_workers=4)

# Structure factor
time_q, q_avg = structure_factor.compute_structure_factor(top, traj, n_frames=n_frames, n_windows=50, num_grid=32, n_workers=4)
```

## Demo
![contact_network_00100](https://github.com/user-attachments/assets/9dd4e0f6-8f2d-4e6e-84b3-6517ea2f62d9)
![contact_network_00200](https://github.com/user-attachments/assets/c45a423c-1787-4a88-8cb1-8b8565231742)
![structure_factor_q](https://github.com/user-attachments/assets/646407be-1d7f-45a3-a9f1-031e2e391573)
![violin_00300](https://github.com/user-attachments/assets/a260b71d-1066-4c0f-bea9-7890315a8ff2)


