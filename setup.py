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
)
