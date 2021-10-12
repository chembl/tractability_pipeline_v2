import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ot_tractability_pipeline_v2",
    version="2.0",
    author="Melanie Schneider, Chris Radoux",
    author_email="melanie@ebi.ac.uk",
    description="Python script for assigning genes to tractability buckets",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/chembl/tractability_pipeline_v2",
    packages=setuptools.find_packages(),
    install_requires=[i.strip() for i in open("requirements.txt").readlines()],
    include_package_data=True,
    entry_points = {
        'console_scripts': ['run-ot-pipeline-v2=ot_tractability_pipeline_v2.bin.run_pipeline:main'],
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
