from setuptools import setup
from setuptools import find_packages

VERSION = '0.1.0'
with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setup(
    name='hellobc',  # package name
    version=VERSION,  # package version
    description='XuLab package for single-cell RNA-seq data upstream processing.',  # package description
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    zip_safe=False,
    project_urls={
        #"Documentation": "https://flask.palletsprojects.com/",
        # "Code": "https://github.com/XintongYao0611/xuscprep",
        #"Issue tracker": "https://github.com/pallets/flask/issues",
    },
    python_requires=">=3.6",
    install_requires=[
        "biopython>=1.80",
        "pysam>=0.18.0",
        "pandas>=1.5.0",
        "numpy>=1.22.0",
        "matplotlib>=3.4.0",
        "scipy>=1.8.0",
        "pyranges>=0.0.85",
    ],
    exclude_package_data={
        "": ["README.md", "LICENSE", "test_pip.py"],
    },
)