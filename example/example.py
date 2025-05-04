if __name__ == '__main__':
    from CALVADOS_Analysis.cluster import size_distribution, time_evolution, largest_size
    from CALVADOS_Analysis.contact import frequency, heatmap, interchain
    from CALVADOS_Analysis.diffusion import compute_msd_and_diffusion
    from CALVADOS_Analysis.structure_factor import compute_structure_factor

    path = r'E:\My_Project\Phase_Seperation\FUS-PLD\CGAnalysis\example\FUS_1_163_500'
    top = f"{path}/checkpoint.pdb"
    traj = f"{path}/FUS_1_163_500_quench_quench.dcd"

    # n_frames = 2000

    # size_distribution.analyze_clusters(top, traj, n_frames=n_frames)
    # time_evolution.analyze_large_clusters(top, traj, n_frames=n_frames)
    # largest_size.compute_largest_size(top, traj, n_frames=n_frames)

    # heatmap.generate_contact_heatmaps(top, traj, n_frames=n_frames)
    # frames, contacts = interchain.compute_interchain_contacts(top, traj, n_frames=n_frames)

    # times, msd, D = compute_msd_and_diffusion(top, traj, n_frames=n_frames, n_workers=20)
    time_q, q_avg = compute_structure_factor(top, traj, n_frames=2000, n_windows=200, num_grid=32, n_workers=15)
