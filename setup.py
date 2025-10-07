#  install eztfem in editable mode: pip install -e .

from setuptools import setup, find_packages

setup(
    name="eztfem",
    version="1.0",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pyvista',
        'ipykernel',
        'ipywidgets',
        'trame',
        'trame-vtk',
        'trame-vuetify'
    ],
    author="Nick Jaensson",
    author_email="n.o.jaensson@tue.nl",
    python_requires='>=3.11',
)
