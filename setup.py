from distutils.core import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'pyreflect',
    packages = ['pyreflect'],
    version = '0.0.1',
    description = "Python codes to work with George Randall's synthetic seismogram code",
    author='Philip Crotwell',
    author_email='crotwell@seis.sc.edu',
    url='https://github.com/crotwell/pyreflect',
    classifiers=["Development Status :: 3 - Alpha",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 3",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
                 "Topic :: Scientific/Engineering",
                 "Operating System :: OS Independent"],
    long_description=long_description,
    long_description_content_type='text/markdown'
)
