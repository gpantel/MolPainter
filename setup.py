import setuptools

setuptools.setup(
    name="MolPainter",
    version="0.9.2",
    description="Tool for drawing complex planar molecular systems of arbitrary composition and molecule placement",
    url='',
    author='George Pantelopulos, Aaron Liberatore',
    license='MIT',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    packages=["MolPainter"],
    package_data={
        "": ['icons/*.png'],
    },
    install_requires=['numpy', 'MDAnalysis'],
    entry_points={
        'console_scripts': [
            'molpainter=MolPainter.painter:main',
        ],
    },
)