import setuptools

setuptools.setup(
    name='libdna',
    version='0.0.4',
    author='Antony B Holmes',
    author_email='antony.b.holmes@gmail.com',
    description='A library for working with DNA.',
    packages=setuptools.find_packages(),
    test_suite='nose.collector',
    tests_require=['nose'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
