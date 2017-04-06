from setuptools import setup

setup(
    name='investor_lifespan_model',
    version='0.0.1',
    author='Nicholas Moehle',
    author_email='nicholasmoehle@gmail.com',
    packages=['investor_lifespan_model'],
    package_dir={'investor_lifespan_model': 'investor_lifespan_model'},
    url='http://github.com/moehle/investor_lifespan_model',
    license='GPLv3',
    description='A simulator of the lifespan of an investor.',
    install_requires=["numpy >= 1.9",
                      "scipy >= 0.15"],
)
