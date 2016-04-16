from setuptools import setup, find_packages

setup(
    name='sequence_tools',
    version='0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        get_seq=sequence_tools.get_seq:get_seq
    	primer_finder=sequence_tools.primer_finder:matching_primer
	    unknown_primer=sequence_tools.unknown_primer:unknown_primer
	    advanced_get_seq=sequence_tools.advanced_get_seq:get_seq
    ''',
)
