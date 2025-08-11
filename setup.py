from distutils.core import setup

setup(
    name='cryspy',
    version='0.0.0',
    author='E. J. Payton',
    author_email='paytonej at ucmail dot uc dot edu',
    packages=['cryspy'],
    scripts=['bin/DeGraefHenry.py',
             'bin/TESTS.py',
             'bin/TUTORIAL.py'],
    url='https://github.com/paytonej/cryspy/',
    license='LICENSE.md',
    description='Visualization and postprocessing of EBSD data.',
    long_description=open('README.md').read(),
    install_requires=[
        'numpy >= 1.0.4',
        'scipy == 0.6.0',
        'platform >= 1.0.7',
        'matplotlib',
        'os',
        'copy',
        'time',
        'sys',
        'h5py'
        ]# 'visvis>= 1.8'
        )