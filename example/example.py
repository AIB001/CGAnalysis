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
    # size_distribution.analyze_clusters(top, traj, n_frames=n_frames)
    # time_evolution.analyze_large_clusters(top, traj, n_frames=n_frames)
    # largest_size.compute_largest_size(top, traj, n_frames=n_frames)

    ##############################
    # Contact Analysis
    ##############################
    # heatmap.generate_contact_heatmaps(top, traj, n_frames=n_frames)
    # frames, contacts = interchain.compute_interchain_contacts(top, traj, n_frames=n_frames)

    ##############################
    # Structural Factor and MSD Analysis
    ##############################
    # time_q, q_avg = compute_structure_factor(top, traj, n_frames=2000, n_windows=200, num_grid=32, n_workers=15)
    # times, msd, D = compute_msd_and_diffusion(top, traj, n_frames=n_frames, n_workers=20)


    #############################
    # Slab Density Calculation
    #############################
    # from CALVADOS_Analysis import density

    # slab_top = os.path.join(base_dir, "slab_case", "checkpoint.pdb")
    # slab_traj =os.path.join(base_dir, "slab_case", "FUS_1_163_500.dcd")

    # z, kde = density.analyze_slab_density(slab_top, slab_traj, bandwidth=0.05, last_n_frames=200, n_workers=10)

