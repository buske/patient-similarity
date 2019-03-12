import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="patient-similarity",
    version="0.0.1",
    author="Orion Buske",
    author_email="buske@cs.toronto.edu",
    description="Phenotype-based patient similarity metric library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/buske/patient-similarity",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'patient-similarity=patient_similarity:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
