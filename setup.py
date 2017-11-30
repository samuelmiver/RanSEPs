# 2017 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'ranseps',
    'description': 'a framework for genome re-annotation and novel small proteins detection adjusting the search to different genomic features that govern protein-coding capabilities',
    'author': 'Samuel Miravet-Verde',
    'url': 'https://github.com/SMV818VMS/RanSEPs',
    'download_url': 'https://github.com/SMV818VMS/RanSEPs',
    'author_email': 'samuel.miravet@crg.edu',
    'version': '0.0.1',
    'license': "GNU GPL",
    'install_requires': ['biopython','numpy','pandas','matplotlib',
                         'scipy','sklearn', 'setuptools'],
    'packages': ['ranseps'],
    'package_data':{'ranseps': ['dbs/*.fa']},
    'scripts': ['bin/ranseps']
}

setup(**config)

# 2017 - Centre de Regulacio Genomica (CRG) - All Rights Reserved
