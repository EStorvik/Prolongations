from setuptools import find_packages, setup

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt") as f:
    required = f.read().splitlines()

with open("requirements-dev.txt") as f:
    required_dev = f.read().splitlines()

setup(
    name="Prolongations",
    version="0.0.1",
    description="For computing symbolic prolongations of vector fields",
    keywords="Lie theory, differential equations, vector fields, prolongations, symmetry",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    package_data={"inverted_pendulum": ["py.typed"]},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python ::3",
        "Programming Language :: Python ::3.10"
        "License :: OSI Approved :: Apache v2 License",
        "Operating System :: Linux, Windows",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=required,
    extras_require={
        "dev": required_dev,
    },
    python_requires=">=3",
    url="https://github.com/EStorvik/Prolongations.git",
    author="Erlend Storvik",
    maintainer="Erlend Storvik",
    maintainer_email="erlend.storvik@hvl.no",
    platforms=["Linux", "Windows"],
    license="GNU General Public License v3",
)
