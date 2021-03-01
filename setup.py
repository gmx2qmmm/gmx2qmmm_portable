from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='gmx2qmmm',
    version='0.0.1',
    packages=setuptools.find_packages(),
    url='https://github.com/gmx2qmmm/gmx2qmmm_portable/tree/p3',
    license='gmx2qmmm',
    author='jangoetze',
    author_email='gmx2qmmm@gmail.com',
    description='A python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=['numpy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
