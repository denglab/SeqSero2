import os, sys
from distutils.core import setup
from setuptools import find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='SeqSero2',
    version=open("version.py").readlines()[-1].split()[-1].strip("\"'"),
    description='Salmonella serotyping',
    long_description=readme(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Programming Language :: Python :: 3',
        'Topic :: Text Processing :: Linguistic',
        ],
    keywords='Salmonella serotyping bioinformatics WGS',
    url='https://github.com/denglab/SeqSero2/',
    author='Shaokang Zhang, Hendrik C Den-Bakker and Xiangyu Deng',
    author_email='zskzsk@uga.edu, Hendrik.DenBakker@uga.edu, xdeng@uga.edu',
    license='GPLv2',
    scripts=["bin/deinterleave_fastq.sh","bin/Initial_Conditions.py","bin/SeqSero2_package.py","bin/SeqSero2_update_kmer_database.py"],
    packages=[""],
    include_package_data = True,
    install_requires=['biopython==1.73'],
    data_files=[("seqsero2_db",["seqsero2_db/antigens.pickle","seqsero2_db/H_and_O_and_specific_genes.fasta","seqsero2_db/invA_mers_dict","seqsero2_db/special.pickle"])],
    zip_safe=False,
)
