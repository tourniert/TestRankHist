import os
from setuptools import setup, find_packages


CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering',
]


setup(name='testrankhist',
      description='Statisctical tests for rank histograms',
      long_description=(open('README.rst').read()
                        if os.path.exists('README.rst')
                        else ''),
      version='0.1',
      license='MIT',
      classifiers=CLASSIFIERS,
      author='tourniert',
      author_email='theo.tournier@meteo.fr',
      url='https://github.com/tourniert/TestRankHist',
      project_urls={
        "Bug Tracker": "https://github.com/tourniert/TestRankHist/issues",
      },
      install_requires=['numpy', 'scipy', 'matplotlib'],
      package_dir={"": "src"},
      packages=setuptools.find_packages(where="src"),
      tests_require=['nose']
)