import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="era2dardar",
    version="0.0.1",
    author="Inderpreet Kaur",
    description="Interpolate ERA5 data to DARDAR resolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SEE-MOF/DARDAR_ERA_interpolation",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU Affero",
        "Operating System :: OS Independent",
    ],
    install_requires=[ "numpy", "xarray", "pyhdf", "typhon", "pansat"],
    python_requires=">=3.6",
    include_package_data=True,
    package_data={
        '': ['*.ini']
    }
)

