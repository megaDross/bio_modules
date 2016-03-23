from setuptools import setup, find_packages

setup(
    name='primer_tools',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        primer_finder=primer_tools.primer_finder:matching_primer
        unknown_primer=primer_tools.unknown_primer:unknown_primer
    ''',
)
