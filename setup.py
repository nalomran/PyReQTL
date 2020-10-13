from setuptools import setup, find_packages

with open('README.md') as f:
    read_me = f.read()

with open('requirements.txt') as rf:
    requirements = rf.read()

setup(
    name='PyReQTL',
    version='0.1.0',
    description='A python library equivalent to R ReQTL Toolkit to '
                'identify the association between expressed SNVs with '
                'their gene expression using RNA-sequencing data.',

    long_description=read_me,
    author='Nawaf Alomran',
    author_email='nawafalomran@hotmail.com',
    url='https://github.com/nalomran/PyReQTL',
    packages=find_packages(),
    install_requires=requirements,
    license='MIT License',
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        'License :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
    python_requires='>=3.5',
)