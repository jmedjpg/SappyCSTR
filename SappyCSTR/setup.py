
from setuptools import setup, find_packages

setup(
    name='SappyCSTR',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pyomo',
        'thermo'
    ],
    author='Jacob Medley',
    author_email='jmedley@andrew.cmu.edu',
    description='A package to optimize CSTRs for saponification reactions.',
    url='https://github.com/jmedjpg/optimize_CSTR',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
)
