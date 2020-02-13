from setuptools import setup
#from setuptools import find_packages

setup(
    name="EpiRank",

    version="0.0.1",

    author="Benny Chin",
    author_email="wcchin.88@gmail.com",

    packages=['EpiRank'],

    include_package_data=False,

    url="https://bitbucket.org/wcchin/epirank3",

    license="LICENSE.txt",
    description="algorithm using forward and backward movements of a commuting network for capturing diffusion.",

    long_description=open("README.md").read(),

    classifiers=[
        'Development Status :: 3 - Alpha',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: GIS',

         'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.6',
    ],

    keywords='spatial network concentration',

    install_requires=[
        "networkx",
        "pandas",
        "numpy",
        "scipy",
    ],
)
