from setuptools import setup, find_packages

with open('README.md') as f:
    read_me = f.read()

with open('requirements.txt') as rf:
    requirements = rf.read()

setup(
    name='PyReQTL',
    version='0.3.0',
    description='A python library equivalent to R ReQTL Toolkit.',
    long_description=read_me,
    long_description_content_type='text/markdown',
    author='Nawaf Alomran',
    author_email='nawafalomran@hotmail.com',
    url='https://github.com/nalomran/PyReQTL',
    download_url='https://github.com/nalomran/PyReQTL/archive/PyReQTL-0.1.0.tar.gz',
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: MIT License",
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
    python_requires='>=3.5',
)