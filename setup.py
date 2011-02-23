from distutils.core import setup

version = __import__('umbrellas').__version__

setup(
    name = "umbrellas",
    version = version.replace(' ', '-'),
    url = 'https://github.com/davecap/umbrellas',
    download_url = 'https://github.com/davecap/umbrellas',
    author = 'David Caplan',
    author_email = 'dcaplan@gmail.com',
    description = 'A small framework for setting up umbrella sampling ensembles',
    packages = ['umbrellas'],
)

