# GlycoGenius: Glycomics Data Analysis Tool
# Copyright (C) 2023 by Hector Franco Loponte
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or 
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. It is accessible within the program files
# or by typing 'license' after running it stand-alone in the terminal
# by typing 'glycogenius'. If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup, find_packages
from distutils.command.install import INSTALL_SCHEMES

for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']
    
long_description_from_file = ""
with open("README.md", "r", encoding="utf-8") as f:
    for lines in f:
        if lines[0] != "!":
            long_description_from_file+= lines
    f.close()

setup(
    name='glycogenius_GUI',
    version='0.0.45',
    author='Hector Franco Loponte',
    author_email='hectorfloponte@gmail.com',
    description='Accessory GUI for Glycogenius',
    long_description=long_description_from_file,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'glycogenius_GUI': ['Assets/*']
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
    install_requires=["pandas==2.1.3", "scipy==1.11.4", "pyteomics==4.6.3",
                     "pillow==10.4.0", "ttkwidgets==0.13.0", "matplotlib==3.8.2",
                     "dill==0.3.7", "numpy==1.26.4", "lxml==4.9.3",
                     "openpyxl==3.1.2", "setuptools==72.2.0",
                     "xlsxwriter==3.2.0", "psutil==5.9.8", "glycogenius==1.1.36"],
    entry_points={
        'gui_scripts': [
            'glycogenius_GUI = glycogenius_GUI:glycogenius_GUI',
        ]
    }
)