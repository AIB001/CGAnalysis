from setuptools import setup, find_packages

setup(
    name="CGAnalysis",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
        "MDAnalysis",
        "mdtraj",
        "pandas",
        "scipy",
        "networkx",
        "seaborn",
        "rdkit",
        "tqdm"
    ],
    author="Zhaoqi SHi",
    description="CGAnalysis: Coarse-grained analysis tools for CALVADOS simulations",
    python_requires=">=3.8"
)
