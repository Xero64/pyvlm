import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pyvlm",
    version="0.0.7",
    author="Xero64",
    author_email="xero64@gmail.com",
    description="Vortex Lattice Method in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Xero64/pyvlm",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    entry_points = {
        'console_scripts': ['pyvlm=pyvlm.__main__:main',],
    }
)
