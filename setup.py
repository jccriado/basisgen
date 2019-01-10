from distutils.core import setup

with open('README.md', mode='r') as readme_file:
    readme = readme_file.read()

setup(
    name='basisgen',
    version='0.3',

    description='Package for the generation of operator bases',
    long_description=readme,

    url='https://github.com/jccriado/basisgen',

    author='Juan Carlos Criado',
    author_email='jccriadoalamo@ugr.es',

    license='MIT',

    packages=['basisgen'],

    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
