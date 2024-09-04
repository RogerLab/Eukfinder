from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="eukfinder",
    version="1.2.3",
    author="dzhao",
    author_email="dandan.tanny.zhao@email.com",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=([
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ]),
    scripts=["bin/Eukfinder.py"],
    entry_points={
        'console_scripts': [
            'eukfinder = Eukfinder:main',
        ],
    }
)
