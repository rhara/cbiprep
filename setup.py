import os, glob
from setuptools import setup, find_packages

setup(
    name='cbiprep',
    version='0.1.9',
    description='CBI Prep',
    author='Ryuichiro Hara',
    author_email='hara.ryuichiro@gmail.com',
    url='https://github.com/rhara/cbiprep',
    license='no license',
    include_package_data=True,
    packages=find_packages(exclude=('tests', 'docs')),
    data_files=[
        ('cbiprep/data', glob.glob('data/*')),
    ],
    entry_points={
        'console_scripts': [
            'ligand_center = cbiprep.ligand_center:main',
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3.7',
    ],
)
