#!/usr/bin/env python
import re
import ast

from pip._internal.req import parse_requirements

from setuptools import setup, find_packages

package_name = 'torsion'

def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]

try:
    install_reqs = get_reqs(parse_requirements("requirements.txt"))
except TypeError:
    from pip._internal.download import PipSession
    install_reqs = get_reqs(
        parse_requirements("requirements.txt", session=PipSession())
    )


setup(
    name=package_name,
    version=__import__(package_name).__version__,
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    author="Vishnu Sresht, Comp Sci",
    author_email="vishnu.sresht@pfizer.com",
    description='Perform Dihedral/Torsional Scans.',
    install_requires=install_reqs,
    license='Other/Proprietary License',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
    ]
)
