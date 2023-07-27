from setuptools import setup

setup(
    name="orphans",
    author="Marina Masson",
    author_email="marina.masson@lpsc.in2p3.fr",
    url = "https://gitlab.in2p3.fr/johan-bregeon/orphans",
    packages=["orphans"],
    description="Studying Gamma-ray Bursts orphan afterglow optical light curves",
    setup_requires=['setuptools_scm'],
    long_description=open("README.md").read(),
    package_data={"": ["README.md", "LICENSE"]},
    use_scm_version={"write_to":"orphans/_version.py"},
    include_package_data=True,
    classifiers=[
        "Development Status :: 1 - Beta",
        "License :: OSI Approved :: GPL License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        ],
    install_requires=["afterglowpy",
                      "matplotlib",
                      "numpy",
                      "astropy",
                      "scipy",
                      "pandas",
                      "pyarrow",
                      "seaborn",
                      "dustmaps",
                      "emcee",
                      "setuptools_scm"]
)
