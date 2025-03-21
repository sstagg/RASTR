from setuptools import setup, find_packages

setup(
    name="RASTR",
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'diameTR=src.scripts.diameTR:main',
			'changestar=src.scripts.changestar:main',
            'averagefft=src.scripts.averagefft:main',
            'RASTR=src.scripts.RASTR_03:main',
            'csexport=src.scripts.csexport:main',
            'scatter_particle_on_micrograph=src.scripts.scatter_particle_on_micrograph:main',
            'showstack=src.scripts.showstack:main',
            'azavg=src.scripts.azavg90:main',
            'createmask=src.scripts.createmask:main',
            'SPOT_RASTR_eliminateduplicates=src.scripts.SPOT_RASTR_eliminateduplicates:main',
            # Add entry point for each script
        ],
    },
    # Add other metadata as needed
    include_package_data=True,
)
