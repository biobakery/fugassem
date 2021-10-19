##### [Replace the contents this README.MD file with the appropriate "Users manual" needed for the python tool.]  
 
# Python Tool Template

## Introduction

> This page gives details concerning guiding principles and formatting required for developing python tool. 

## Requirements
- #### Git and Github accounts
    Follow this installation guide [HERE](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) to install Git in your local machine. Additionally, direct to Github [sign up](https://github.com/join?source=header-home) page if you do not have the account. 
- #### PIP
    Follow [official PIP link](https://pip.pypa.io/en/stable/installing/) for installation instruction. 

     
## Installation

> Getting started with the Python Tool Template:  
- Click on "**Use this Template**" button and fill out the new repository informations (name/description/access). 
- Finally click on "**Create repository from template**" and you now have a new repository based on the Python tools Template following the standard layouts. 

Alternatively, 
- Direct to [Github Create New Repository](https://github.com/organizations/biobakery/repositories/new) and select the "**biobakery/python-tools-template**". Fill out other information (Repository name/description/access)
- Finally click on "**Create repository**" and you now have a new repository based on the Bioconductor Template. 

**NOTE**: Make your repository "**Private**" unless it is ready to be released.

## Getting started with a public python tools

- Clone your recently created repository in your local development environment either using:
    ``` 
        git clone https://github.com/biobakery/Python-tool-template.git
    ```
    or using the "**Clone or Download**" button. 

### Picking A Name
Python module/package names should generally follow the following constraints:
- All lowercase
- Unique on pypi, even if you don’t want to make your package publicly available (you might want to specify it privately as a dependency later)
- Underscore-separated or no word separators at all (don’t use hyphens)

### Packaging
The main setup config file, **setup.py**, should contain a single call to setuptools.setup(), like so:
```
from setuptools import setup
setup(name='<Package Name>',
      version='<Version>',
      description='<Description of the package>',
      author='Author Name',
      author_email='Author Email',
      license='MIT',
      packages=['<package name>'],
      zip_safe=False)
```
Addition parameter included in the sample **setup.py** template file. 
Now we can install the package locally (for use on our system), with:
~~~
$ pip install .
~~~
We can also install the package with a symlink, so that changes to the source files will be immediately available to other users of the package on our system:
~~~
$ pip install -e .
~~~
Anywhere else in our system using the same Python, we can do this now:
~~~
>>> from template.tools import sample_function
>>> print sample_function.transpose()
~~~
**Note:** Edit the **setup.py** file with relevent information about the package. 

For more packaging information, see [official python packaging documentation](https://python-packaging.readthedocs.io/en/latest/)

### Test workflow 
 Two of the main frameworks for testing python tools are [Unittests](https://docs.python.org/3/library/unittest.html) and [pytest](https://docs.pytest.org/en/latest/). 
 
 For this template, we are using unittest: 
 - To set up your package to use unittest, run:
    Either use the module and class name
    ```
    python -m unittest template.tests.sample_test.TestStringMethods
    ```
    Or use the path: 
    ```
    python template/tests/sample_test.py 
    ```
    
    This will run your test cases of the class or the test file. For more information on running your test cases through command line, see [Python unit test documentation](https://docs.python.org/3/library/unittest.html).
