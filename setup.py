from setuptools import setup, find_packages

setup(
    name="bands_NMREGIS",
    version="0.0.1",
    py_modules=["bands_NMREGIS"],
    author="Natan Moreira Regis",
    description="Python package for plot VASP's bands structure",
    classifiers=["Programming Language :: Python :: 3"],
    packages=find_packages(),
    install_requires=[
        "matplotlib",
        "pandas", 
        "numpy"
    ],
    entry_points={
        'console_scripts': [
            'bands_NMREGIS = bands_NMREGIS:menu',  # assumes a `main()` function
        ],
    },
)