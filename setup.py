from setuptools import setup, find_packages


# Long description
with open("README.md", "r") as readme:
    desc = readme.read()

# Short description
sdesc = ""

# Dependencies
with open("requirements.txt") as f_:
    requirements = f_.read().strip().split("\n")

setup(
    name="qmx2qmmm",
    version="1.0.2",
    package_dir={
        'gmx2qmmm': 'gmx2qmmm',
        'gmx2qmmm.operations': 'gmx2qmmm/operations',
        'gmx2qmmm.pointcharges': 'gmx2qmmm/pointcharges'
        },
    packages=[
        'gmx2qmmm',
        'gmx2qmmm.operations',
        'gmx2qmmm.pointcharges'
        ],
    # packages=find_packages(),
    author="Jan Goetze",
    author_email="",
    description=sdesc,
    long_description=desc,
    long_description_content_type="text/markdown",
    url="https://github.com/gmx2qmmm/gmx2qmmm_portable",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ]
)
