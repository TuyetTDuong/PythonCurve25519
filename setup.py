from __future__ import print_function
from glob import glob
from setuptools import setup,Extension

sources = ['source/curve25519module.c', 'source/curve25519_dh.c']
sources.extend(glob("source/curve25519_mehdi.c"))
sources.extend(glob("source/curve25519_order.c"))
sources.extend(glob("source/curve25519_utils.c"))
sources.extend(glob("source/ed25519_sign.c"))
sources.extend(glob("source/ed25519_verify.c"))
sources.extend(glob("source/sha512.c"))
sources.extend(glob("source/custom_blind.c"))
module_curve = Extension('curve25519_topl',
                         sources = sorted(sources),
                         #                   headers = headers,
                         include_dirs = [
                             'source'
                         ]
                         )
setup(
    name='curve25519_topl',
    version="0.4.1-2",
    license='GPLv3 License',
    author='Tuyet Duong',
    ext_modules = [module_curve],
    author_email='duongtt3@alumni.vcu.edu',
    description='curve25519',
    platforms='any'
)