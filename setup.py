import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='HLAcc',
    version='1.0.0',
    author='Yue (Shawn) Shen',
    author_email='yshen25@vols.utk.edu',
    description='Package for HLA clustering based on coarse-grained structures',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yshen25/HLAcc',
    project_urls={
        'Documentation': 'https://github.com/yshen25/HLAcc',
        'Bug Reports':
        'https://github.com/yshen25/HLAcc/issues',
        'Source Code': 'https://github.com/yshen25/HLAcc'
    },
    package_dir={'': 'src'},
    # packages=['HLAcc'],
    packages=setuptools.find_packages(where='src'),
    # package_data={'HLAcc': ['dat/*']},
    classifiers=[
        # see https://pypi.org/classifiers/
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Healthcare Industry',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent'
    ],
    python_requires='>=3.8, <3.9',
    install_requires=['numpy', 'scipy', 'pandas', 'biopandas', 'biopython', 'matplotlib', 'seaborn', 'scikit-learn', 'pymol']
)
