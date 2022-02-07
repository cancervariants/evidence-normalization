"""Module for package and distribution"""
from setuptools import setup

exec(open("evidence/version.py").read())
setup(version=__version__)
