import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eukfindertest2",
    version="1.0.1",
    author="dzhao",
    author_email="dandan.tanny.zhao@email.com",
    description="A small demo package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dzhao2019/eukfindertest",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
