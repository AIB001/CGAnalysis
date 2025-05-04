import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fftn, fftfreq, fftshift
from tqdm import tqdm
from scipy.stats import linregress
import multiprocessing as mp
import os

def _process_window(args):
    win, top, traj, num_grid, window_size, n_windows = args

    u = mda.Universe(top, traj)
    ca = u.select_atoms('name CA')

    start = win * window_size
    end = min((win + 1) * window_size, len(u.trajectory))

    box = u.trajectory.ts.dimensions

    rho_accum = np.zeros((num_grid, num_grid, num_grid), dtype=np.float32)

    for ts in u.trajectory[start:end]:
        pos = ca.positions - ts.dimensions[:3] / 2
        hist, _ = np.histogramdd(pos, bins=num_grid, range=[[-box[0]/2, box[0]/2]]*3)
        rho_k = fftn(hist)
        power = np.abs(rho_k)**2 / len(ca)
        rho_accum += power

    S_q = rho_accum / window_size

    qx = fftfreq(num_grid, d=box[0]/num_grid) * 2*np.pi
    qy = fftfreq(num_grid, d=box[1]/num_grid) * 2*np.pi
    qz = fftfreq(num_grid, d=box[2]/num_grid) * 2*np.pi

    q_mags = np.sqrt(qx[:, None, None]**2 + qy[None, :, None]**2 + qz[None, None, :]**2)
    q_mags = fftshift(q_mags)
    S_q_shifted = fftshift(S_q)

    numerator = np.sum(q_mags * S_q_shifted)
    denominator = np.sum(S_q_shifted)
    q_avg = numerator / denominator if denominator > 0 else 0.0

    return q_avg

def compute_structure_factor(top, traj, n_frames=2000, n_windows=200, num_grid=32, n_workers=4):
    """
    Compute ⟨q⟩ over time using structure factor analysis with parallel processing.
    """

    u = mda.Universe(top, traj)
    total_frames = min(n_frames, len(u.trajectory))
    window_size = total_frames // n_windows

    frame_time_ps = u.trajectory.dt  # assume time unit is ps

    print(f"Using {n_workers} CPU cores for structure factor calculation")

    args_list = [(win, top, traj, num_grid, window_size, n_windows) for win in range(n_windows)]

    with mp.Pool(processes=n_workers) as pool:
        q_avg_list = list(tqdm(pool.imap(_process_window, args_list),
                               total=n_windows))

    q_avg_list = np.array(q_avg_list)
    time = np.arange(n_windows) * window_size * frame_time_ps

    ## ---- Plotting ----
    import seaborn as sns
    sns.set_style("whitegrid")
    plt.rcParams.update({'font.size': 20, 'font.family': 'Times New Roman'})

    plt.figure(figsize=(10, 7))
    plt.plot(time, q_avg_list, 'o-', color='#4370B4', label='Data')

    ## ---- Best fit ----
    n_points = len(q_avg_list)
    min_window = 50
    best_slope = 0
    best_fit = None

    for start in range(0, n_points - min_window, 5):
        for end in range(start + min_window, n_points, 5):
            log_time = np.log(time[start:end] + 1e-6)
            log_q = np.log(q_avg_list[start:end] + 1e-6)
            slope, intercept, r_value, _, _ = linregress(log_time, log_q)
            if abs(slope) > abs(best_slope) and r_value**2 > 0.7:
                best_slope = slope
                best_fit = {
                    'start': start, 'end': end,
                    'slope': slope, 'intercept': intercept,
                    'r_value': r_value
                }

    if best_fit:
        start, end = best_fit['start'], best_fit['end']
        fit_time = time[start:end]
        fit_line = np.exp(best_fit['intercept']) * fit_time**best_fit['slope']

        extension = 0.1 * (fit_time[-1] - fit_time[0])
        extended_time = np.linspace(max(fit_time[0] - extension * 8, 1), fit_time[-1] + extension, 100)
        extended_fit_line = np.exp(best_fit['intercept']) * extended_time**best_fit['slope']

        plt.plot(extended_time, extended_fit_line, 'r--', linewidth=3,
                 label=f'Best fit: slope={best_fit["slope"]:.3f} (R²={best_fit["r_value"]**2:.3f})')
        plt.plot(fit_time, fit_line, 'r-', linewidth=3)
        plt.axvspan(time[start], time[end-1], color='#549F9A', alpha=0.2)

        ## ---- Fixed slope -1/3 ----
        log_time_fixed = np.log(fit_time)
        log_q_fixed = np.log(q_avg_list[start:end] + 1e-6)
        fixed_slope = -1/3
        fixed_intercept = np.mean(log_q_fixed) - fixed_slope * np.mean(log_time_fixed)

        fixed_fit_line = np.exp(fixed_intercept) * fit_time**fixed_slope
        fixed_extended_fit_line = np.exp(fixed_intercept) * extended_time**fixed_slope

        residuals = log_q_fixed - (fixed_intercept + fixed_slope * log_time_fixed)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((log_q_fixed - np.mean(log_q_fixed))**2)
        r_squared = 1 - (ss_res / ss_tot)

        plt.plot(extended_time, fixed_extended_fit_line, 'g--', linewidth=3,
                 label=f'Fixed slope=-1/3 (R²={r_squared:.3f})')
        plt.plot(fit_time, fixed_fit_line, 'g-', linewidth=3)

    plt.xlabel('Time (ps)', fontname='Times New Roman')
    plt.ylabel(r'$\langle q \rangle$ ($\AA^{-1}$)', fontname='Times New Roman')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(alpha=0.3)
    plt.legend()

    out_dir = os.path.dirname(traj)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "structure_factor_q.png"), dpi=300)
    plt.close()

    return time, q_avg_list
