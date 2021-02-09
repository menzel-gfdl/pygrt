from glob import glob
from os.path import join
from setuptools import Extension, find_packages, setup


def grtcode():
    dirs = ["GRTCODE/{}/src".format(x) for x in ["utilities", "gas-optics"]]
    src = []
    for directory in dirs:
        src += glob(join(directory, "*.c"))
    return Extension("pygrt.libgrtoptics",
                     sources=src,
                     include_dirs=dirs,
                     extra_compile_args=["-fopenmp"],
                     extra_link_args=["-fopenmp"])


setup(
    name="pygrt",
    version="0.0.0",
    author="R. Menzel",
    author_email="author@example.com",
    description="Python frontend for GRTCODE.",
    url="",
    python_requires=">=3.5",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: LGPL-2.1 License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy", 
        "pyrad @ git+http://github.com/menzel-gfdl/pylbl@line-mixing",
    ],
    ext_modules = [grtcode()],
)
