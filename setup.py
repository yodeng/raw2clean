#!/usr/bin/env python
from setuptools import setup

setup(
    name = "raw2clean",
    version = "1.0.8",
    packages = ["raw2clean"],
    package_dir = {"raw2clean":"src"},
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

