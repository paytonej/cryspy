from distutils.core import setup

setup(
    name='cryspy',
    version='0.0.0',
    author='E. J. Payton',
    author_email='eric dot payton at us dot af dot mil',
    packages=['cryspy'],
    scripts=['bin/DeGraefHenry.py',
             'bin/TESTS.py',
             'bin/TUTORIAL.py'],
    url='http://www.wpafb.af.mil/afrl/rx/',
    license='LICENSE.txt',
    description='Visualization and postprocessing of EBSD data.',
    long_description=open('README.txt').read(),
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