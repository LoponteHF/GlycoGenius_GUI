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
import sys
import os

current_version='1.0.3'

for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

def bump_version(version, bump_type):
    current_version = version
    major, minor, patch = map(int, current_version.split("."))

    if bump_type == "major":
        major += 1
        minor = 0
        patch = 0
    elif bump_type == "minor":
        minor += 1
        patch = 0
    elif bump_type == "patch":
        patch += 1
    else:
        return version

    new_version = f"{major}.{minor}.{patch}"
    return new_version

# Bumps version
if len(sys.argv) > 1 and 'sdist' not in sys.argv:
    bump_type = sys.argv[1]
    if bump_type in ["major", "minor", "patch"]:
        new_version = bump_version(current_version, bump_type)
        with open(__file__, "r") as file:
            lines = file.readlines()
        with open(__file__, "w") as file:
            for line in lines:
                if line.startswith("current_version="):
                    file.write(f"current_version='{new_version}'\n")
                else:
                    file.write(line)
        print("Version bumped succesfully!")
        print(f"Former version: {current_version}, New version: {new_version}")
        os._exit(0)

# Captures the description from the markdown readme    
long_description_from_file = ""
with open("README.md", "r", encoding="utf-8") as f:
    for lines in f:
        if lines[0] != "!":
            long_description_from_file+= lines
    f.close()

# Captures the requirements from the requirements file
requirements = []
with open("requirements.txt", "r") as f:
    for line in f:
        if line[0] != "#":
            requirements.append(line.strip())

setup(
    name='glycogenius_GUI',
    version=current_version,
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
    install_requires=requirements,
    entry_points={
        'gui_scripts': [
            'glycogenius_GUI = glycogenius_GUI:glycogenius_GUI',
        ]
    }
)