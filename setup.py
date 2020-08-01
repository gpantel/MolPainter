import setuptools

setuptools.setup(
    name="MolPainter",
    version="0.9.3",
    description="Tool for drawing complex planar molecular systems of arbitrary composition and molecule placement",
    url='https://github.com/gpantel/MolPainter',
    author='George Pantelopulos, Aaron Liberatore',
    license='MIT',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    packages=["MolPainter", "MolSolvator"],
    package_data={
        "": ['icons/*.png'],
    },
    install_requires=['numpy', 'scipy', 'toml', 'MDAnalysis'],
    entry_points={
        'console_scripts': [
            'molpainter=MolPainter.painter:main',
            'molsolvator=MolSolvator.MolSolvator:main',
        ],
    },
)