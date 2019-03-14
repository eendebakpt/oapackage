# Contributing 

### Code style

* Try to follow the [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guide. A usefull tool for automated formatting is [autopep8](https://pypi.python.org/pypi/autopep8). We do allow lines upto 120 characters.
* Document functions using the [Google docstring](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html) convention
* Add unit tests using the [unittest](https://docs.python.org/3/library/unittest.html) framework

### Submitting code

If you would like to contribute, please submit a pull request. (See the [Github Hello World](https://guides.github.com/activities/hello-world/) example, if you are new to Github).

By contributing to the repository you state you own the copyright to those contributions and agree to include your contributions as part of this project under the BSD license.

### Bugs reports and feature requests

To submit a bug report or feature request use the [Github issue tracker](https://github.com/eendebakpt/oapackage/issues).
Search for existing and closed issues first. If your problem or idea is not yet addressed, please open a new issue.

### Testing

Continuous integration and testing for the C++ library is performed on [Travis](https://travis-ci.org/eendebakpt/oapackage) and for the Python package on [AppVeyor](https://ci.appveyor.com/project/eendebakpt/oapackage-4lws8).

To perform tests run [`pytest`](https://docs.pytest.org/en/latest/). To obtain a [coverage](https://coverage.readthedocs.io) report, run
```
$ coverage run --source='./oapackage' -m pytest
$ coverage report --omit oapackage/markup.py,oapackage/_scanf.py
```

### Contact

For further information please contact Pieter Eendebak, pieter.eendebak at tno.nl

