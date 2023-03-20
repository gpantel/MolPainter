import setuptools
import MolPainter

setuptools.setup(
    name='MolPainter',
    version=MolPainter.__version__,
    description='Tool for drawing complex planar molecular systems of arbitrary composition and molecule placement',
    long_description=open('README.md', 'r').read(),
    long_description_content_type='text/markdown',
    url=MolPainter.__url__,
    author='George Pantelopulos, Aaron Liberatore',
    license='MIT',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    packages=['MolPainter'],
    package_data={
        '': ['icons/*.png'],
    },
    install_requires=['numpy', 'scipy', 'toml', 'MDAnalysis'],
    entry_points={
        'console_scripts': [
            'molpainter=MolPainter.painter:main',
            'molsolvator=MolPainter.MolSolvator:main',
        ],
    },
)
