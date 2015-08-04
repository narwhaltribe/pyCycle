from setuptools import setup

kwargs = {'author': '',
 'author_email': '',
 'classifiers': ['Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering'],
 'data_files': [('pycycle', ['src/pycycle/gri1000.cti'])],
 'description': '',
 'download_url': '',
 'include_package_data': True,
 'install_requires': [],
 'keywords': ['openmdao'],
 'license': '',
 'maintainer': '',
 'maintainer_email': '',
 'name': 'pycycle',
 'package_data': {'pycycle': []},
 'package_dir': {'': 'src'},
 'packages': ['pycycle'],
 'url': '',
 'version': '0.1',
 'zip_safe': False}


setup(**kwargs)

