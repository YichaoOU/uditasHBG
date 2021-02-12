# -*- coding: utf-8 -*-


"""setup.py: setuptools control."""

# Followed the advice in https://github.com/jgehrcke/python-cmdline-bootstrap

from setuptools import setup

from uditasHBG._version import __version__

version = __version__

with open("README.md", "rb") as f:
    long_descr = f.read().decode("utf-8")


setup(
    name = "uditasHBG",
    packages = ["uditasHBG"],
    entry_points = {
        "console_scripts": ['uditasHBG = uditasHBG.uditas:main']
        },
    version = version,
    description = "UDiTaS analysis software.",
    long_description = long_descr,
    author = "Editas Medicine, Inc",
    url = "http://www.editasmedicine.com/",
    zip_safe = False,
    )
