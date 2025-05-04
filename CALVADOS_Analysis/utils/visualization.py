import matplotlib.pyplot as plt
import seaborn as sns

def save_plot(x, y, xlabel, ylabel, title, filename, line_style='-o'):
    """
    Save a simple line plot.

    Parameters:
    - x: X-axis data.
    - y: Y-axis data.
    - xlabel: X-axis label.
    - ylabel: Y-axis label.
    - title: Plot title.
    - filename: Output filename.
    - line_style: Line style for the plot.
    """
    sns.set_style("whitegrid")
    plt.figure()
    plt.plot(x, y, line_style)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
