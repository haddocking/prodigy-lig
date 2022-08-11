import os
from setuptools import setup
import codecs

from prodigy_lig import version


def read(fname):
    return codecs.open(os.path.join(os.path.dirname(__file__), fname), encoding='utf-8').read()

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name='prodigy_lig',
    version = version,
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
    install_requires=required,
    entry_points={
        'console_scripts': [
            'prodigy_lig = prodigy_lig.prodigy_lig:main',
        ]
    },
    zip_safe=False
)
