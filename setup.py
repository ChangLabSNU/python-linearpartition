#
# Copyright 2023 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# “Software”), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

from setuptools import setup, Extension
import numpy
import os

numpy_includes = [numpy.get_include()]

lpartitionext = Extension(
            'linearpartition',
            ['linearpartitionmodule.cc'],
            include_dirs=['LinearPartition/src', 'contrib'] + numpy_includes)

if not os.path.exists('LinearPartition/src/LinearPartition.cpp'):
    print('''
## The source code of LinearPartition is required to build this module.
## Please download the code by running the following command in your terminal:
##   git clone https://github.com/LinearPartition/LinearPartition.git
''')

setup(
    name='linearpartition-unofficial',
    version='0.2',
    description='Python interface to LinearPartition, a linear-time RNA secondary structure prediction tool',
    author='Hyeshik Chang',
    author_email='hyeshik@snu.ac.kr',
    url='https://github.com/ChangLabSNU/python-linearpartition',
    download_url='https://github.com/ChangLabSNU/python-linearpartition/releases',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    keywords=[
        'RNA',
        'secondary structure',
        'RNA structure'
    ],
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Plugins',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    ext_modules=[lpartitionext],
    setup_requires=['numpy'],
    install_requires=['numpy'],
)
