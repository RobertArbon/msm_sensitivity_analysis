import setuptools

setuptools.setup(
    name="sensetools",
    version="0.1",
    description="MSM sensitivity analysis",
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "addict",
        "click",
        "devtools",
        "mdtraj",
        "numpy",
        "deeptime",
        "h5py",
        "pretty-errors",
        "pandas",
        "tables",
        "matplotlib"
    ],
    entry_points={
        'console_scripts': [
            'sensetools = sensetools.main:cli',
        ],
    },
    python_requires=">=3.8",
)