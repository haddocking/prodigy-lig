import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

requirements = [
    "numpy",
    "biopython",
]

setup(
    name='prodigy_lig',
    version = "0.0.0",
    description=("Calculate protein-small molecule binding affinities."),
    url='http://github.com/haddocking/prodigy-lig',
    author='Computational Structural Biology Group @ Utrecht University',
    author_email='prodigy.bonvinlab@gmail.com',
    license='Apache 2.0',
    packages=['prodigy_lig'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 4 - Beta",
    ],
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'prodigy_lig = prodigy_lig.prodigy_lig:main',
        ]
    },
    zip_safe=False
)
