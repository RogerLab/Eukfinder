from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="eukfinder",
    version="1.2.4",
    author="Dayana E. Salas-Leiva, Dandan Zhao",
    author_email="ds2000@cam.ac.uk, d.zhao@dal.ca",
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
