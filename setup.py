from setuptools import setup

setup(
    name='aln_stats',
    version='0.1',
    py_modules=['aln_stats'],
    entry_points={
        'console_scripts': [
            'aln_stats = aln_stats:main',
        ],
    },
)