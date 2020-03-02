import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(

     name='sasa_phys',

     version='0.1',

     scripts=[] ,

     author="Tim Luca Turan, Max Bräuer",

     author_email="timturan@web.de",

     description="Implementation of the Semi Analytic Stacking Algorithm for Meta Surface stacks",

     long_description=long_description,

   long_description_content_type="text/markdown",

     url="https://github.com/TimLucaTuran/SASA",

     packages=setuptools.find_packages(),

     classifiers=[

         "Programming Language :: Python :: 3",

         "License :: OSI Approved :: MIT License",

         "Operating System :: OS Independent",

     ],

 )
