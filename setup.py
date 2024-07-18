from setuptools import setup

setup(
    name="pyprism",
    version="0.0.1",
    description="""A deconvolution framework for bulk RNA-seq count data based on BayesPrism, implementing 
    the improved algorithm of InstaPrism in Python, integrated into the ScVerse, and interoperable with AnnData and 
    HDF5 files.""",
    url="https://github.com/LuJoHae/PyPrism",
    author="Lukas Jonathan Haeuser",
    author_email="lukas.haeuser@empa.ch",
    license='MIT',
    packages=['pyprism'],
    install_requires=["numpy==1.26",
                      "beartype==0.18.5",
                      "anndata==0.10.8",
                      "h5py==3.11.0",
                      "pyensembl==2.3.13",
                      "typing_extensions==4.12.2"
                      ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.12',
    ],
)