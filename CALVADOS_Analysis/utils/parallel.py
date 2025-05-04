import multiprocessing as mp

def parallel_map(func, iterable, processes=None):
    """
    Apply a function in parallel over an iterable.

    Parameters:
    - func: Function to apply.
    - iterable: Iterable of input values.
    - processes: Number of parallel processes.

    Returns:
    - List of results.
    """
    if processes is None:
        processes = max(1, mp.cpu_count() - 1)
    with mp.Pool(processes) as pool:
        results = pool.map(func, iterable)
    return results
