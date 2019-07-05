import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

data_files = ['README.md', 'LICENSE', 'documentation/Equations.pdf',
              '*.ipynb', 'validation_data/*.txt', 'wptherml/datalib/*.txt']

setuptools.setup(
    include_package_data=True,
    name="wptherml",
    #packages=['wptherml'],
    packages=setuptools.find_packages(),
    package_data={'wptherml':data_files},
    package_dir={'wptherml': 'wptherml'},
    #packages=setuptools.find_packages(),
    #package_dir={'wptherml': 'wptherml'},
    #package_data={'wptherml': ['datalib/*.txt', 'validation_data/*.txt']},
    version="1.0.12b",
    author="Foley Lab",
    author_email="foleyj10@wpunj.edu",
    description="A Python package for the design of materials for harnessing heat.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://foleylab.github.io/wptherml/",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: MacOS :: MacOS X",
	"Operating System :: Microsoft :: Windows :: Windows 10",
	"Operating System :: POSIX :: Linux",
    ],
)
