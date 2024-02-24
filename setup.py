# --------------------------------------------------------------------------#
#								pySPHinXsys									#
# --------------------------------------------------------------------------#
# pySPHinXsys provides a python binding for SPHinXsys library.				#
#																			#
# Copyright(c) 2017-2024 The authors and the authors' affiliations.         #
#                                                                           #
# Licensed under the Apache License, Version 2.0 (the "License"); you may   #
# not use this file except in compliance with the License. You may obtain a #
# copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        #
#                                                                           #
# --------------------------------------------------------------------------#
# @file 	setup.py
# @brief 	Setup for pySPHinXsys stallation. 
# @author	Chi ZHANG
# @version  1.0 
#			Try to implement pybind for SPHinXsys library. 
#			-- Chi ZHANG

import os
import re
import sys
import platform
import subprocess
import multiprocessing as mp
import argparse

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

# Extract cmake arguments
parser = argparse.ArgumentParser()
parser.add_argument("-D", action='append', dest='cmake',
                    help="CMake Options")
args, other_args = parser.parse_known_args(sys.argv)
cmake_clargs = args.cmake
sys.argv = other_args

# Project binding name
name = "pysphinxsys"


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.16.0':
                raise RuntimeError("CMake >= 3.16.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        bin_dir_windows = os.path.join(os.path.abspath(self.build_temp), "bin")
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Release'
        build_args = ['--config', cfg]

        # Add cmake command line arguments
        if cmake_clargs is not None:
            cmake_args += ['-D{}'.format(arg) for arg in cmake_clargs]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                           '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + bin_dir_windows,
                           '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), bin_dir_windows)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j{}'.format(mp.cpu_count())]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        # Add position independent code flags if using gcc on linux probably
        if platform.system() == "Linux":
            cmake_args += ['-DCMAKE_C_COMPILER_LAUNCHER=ccache',
                           '-DCMAKE_CXX_COMPILER_LAUNCHER=ccache', 
                           '-DSPHINXSYS_USE_FLOAT=OFF',
                           '-DSPHINXSYS_MODULE_OPENCASCADE=ON',
                           '-DCMAKE_INSTALL_RPATH={}'.format("$ORIGIN"),
                           '-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON',
                           '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=OFF']

        if platform.system() == "Darwin":
            cmake_args += ['-DCMAKE_C_COMPILER_LAUNCHER=ccache',
                           '-DCMAKE_CXX_COMPILER_LAUNCHER=ccache', 
                           '-DCMAKE_INSTALL_RPATH={}'.format("$ORIGIN"),
                           '-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON',
                           '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=OFF']

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.', '--target', name] + build_args, cwd=self.build_temp)

setup(
    name=name,
    version='1.0.0',
    author='Chi ZHANG',
    author_email='zhangchi0118@gmail.com',
    description='Python Bindings for SPHinXsys library',
    long_description='SPH Simulations using the SPHinXsys package natively from python',
    ext_modules=[CMakeExtension(name)],
    cmdclass=dict(build_ext=CMakeBuild),
    packages=find_packages(),
    zip_safe=False,
    python_requires=">=3.16",
)

#python3 setup.py bdist_wheel -DCMAKE_TOOLCHAIN_FILE="/home/xyhu/vcpkg/scripts/buildsystems/vcpkg.cmake"