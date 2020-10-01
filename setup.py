import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Exo-DMC-mbonav", # Replace with your own username
    version="X.YaN",
    author="Mariangela Bonavita",
    author_email="mbonav@roe.ac.uk",
    description="Exoplanet Detection Map Calculator",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mbonav/Exo-DMC",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
