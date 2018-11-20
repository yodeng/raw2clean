#!/usr/bin/env python
from setuptools import setup, find_packages  

setup(
    name = "raw2clean",
    version = "1.0.7",
    packages = find_packages(),
    author="Yong Deng",
    author_email = "yodeng@tju.edu.cn",
    description = "QC for rawdata without adapter",
    url="https://github.com/yodeng/raw2clean",
    license="MIT",
    entry_points = {
        'console_scripts': [  
            'raw2clean = raw2clean.raw2clean:main',
            'sumfq = raw2clean.summary_fastq:main'
        ]
    }

)

