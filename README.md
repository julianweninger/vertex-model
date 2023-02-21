# Utopia Model Repository

This private repository contains the vertex models within [Utopia] developed in
Karsten Kruse's group at University of Geneva and Madan Rao's group at Simon
Centre, National Centre for Biological Sciences, Bangalore.

**Note:** If not mentioned explicitly, all instructions and considerations from the [main repository][Utopia] still apply here. ☝️

#### Contents of this Readme
* [Installation](#installation)
* [Model Documentation](#model-documentation)
* [Quickstart](#quickstart)
* [Information for Developers](#information-for-developers)
* [Dependencies](#dependencies)

---

## Installation
The following instructions assume that you built Utopia in a
development environment, as indicated in the framework repository's
[`README.md`](https://gitlab.com/utopia-project/utopia/-/blob/master/README.md)

### Step-by-step Instructions
These instructions are intended for 'clean' __Ubuntu__ or __macOS__ setups
with Utopia already _built_ by your user.

#### 0 — Setup Utopia
Make sure to `git checkout` your desired branch of the main Utopia repository,
and to properly configure and build the package as described in its README.

#### 1 — Cloning This Repository
Enter the `Utopia` directory into which you previously cloned the main
repository.

With access to the Utopia GitLab group, you can clone the repository to that
directory using the following command:

```bash
git clone https://gitlab.com/utopia-project/utopia.git
```

Inside your top level `Utopia` directory, there should now be two repositories:
* `utopia`: the framework code
* `vertex-model`: *this* repository, with the vertex model's code and associated
                  model's code

#### 2 — Install dependencies
Install the third-party dependencies using a package manager. We recommend
APT on Ubuntu and Homebrew on macOS. Please consider the notes given in the
section on installing dependencies in the Utopia framework repository.

No additional dependencies required.

#### 3 — Configure and build
Enter the repository and create your desired build directory:

```bash
cd vertex-model
mkdir build
```

Now, enter the build directory and invoke CMake:

```bash
cd build
cmake ..
```

The terminal output will show the configuration steps, which includes the
installation of further Python dependencies and the creation of a virtual
environment.

After this, you can build a specific or all vertex-models using:

```bash
make PCPVertex   # builds only the PCPVertex model
make -j4 all   # builds all models, using 4 CPUs
```

#### 4 — Run a model :tada:
You should now be able to run Utopia models from both the framework and this
repository via the Utopia CLI.

Upon configuration, CMake creates symlinks to the Python virtual environment
of the Utopia repository. You can enter it by navigating to either of the build
directories and calling

```bash
source ./activate
```

You can now execute the Utopia frontend as you are used to, only now the
models of this repository should be available. List all available models
with

```bash
utopia models ls
```

and perform a simulation run for e.g. `PCPTopology` with

```bash
utopia run PCPTopology --cfg-set proliferation
```

#### 5 — Stay up-to-date
This last step remains relevant even after you finished your installation.
Have a look at the [contribution guide](CONTRIBUTING.md) to see the tasks that you should carry out periodically in order to stay up-to-date with both your installation and your model implementation.


### Troubleshooting
* If the `cmake ..` command fails because it cannot locate the Utopia main
    repository, something went wrong with the CMake User Package Registry.
    You can always specify the location of the package _build_ directory via
    the CMake variable `Utopia_DIR` and circumvent this CMake feature:

    ```bash
    cmake -DUtopia_DIR=<path/to/>utopia/build ..
    ```

* If your model requires a certain branch of the Utopia framework, make sure
    to re-configure and re-build the framework repository after switching
    the branch.



## Model Documentation
Unlike the [main Utopia documentation](https://docs.utopia-project.org/html/index.html), the models included in this repository come with their own documentation which has to be built locally.
It is *not* available online.

To build these docs locally, navigate to the `build` directory and execute

```bash
make doc
```

The [Sphinx](http://www.sphinx-doc.org/en/master/index.html)-built user
documentation will then be located at `build/doc/html/`, and the
C++ [doxygen](http://www.stack.nl/~dimitri/doxygen/)-documentation can be
found at `build/doc/doxygen/html/`. Open the respective `index.html` files to
browse the documentation.



## Quickstart
### How to run a model?
The Utopia command line interface (CLI) is, by default, only available in a
Python virtual environment, in which `utopya` (the Utopia frontend) and its
dependencies are installed. This virtual environment is located in the
**main** (`utopia`) repository's build directory. However, symlinks to the
`activate` and `run-in-utopia-env` scripts are provided within the build
directory of this project. You can enter it by specifying the correct path:

```bash
source <path/to/{utopia,vertex-model}>/build/activate
```

Now, your shell should be prefixed with `(utopia-env)`. All the following
should take place inside this virtual environment.

As you have already done with the `dummy` model, the basic command to run a
model named `SomeModel` is:

```bash
utopia run SomeModel
```

You can list all models registered in the frontend with

```bash
utopia model ls
```



## Information for Developers
### New to Utopia? How do I implement a model?
Please refer to the [documentation of the main repository](https://hermes.iup.uni-heidelberg.de/utopia_doc/latest/html/) :wink:

Note that your *personal* model should be implemented in *this* repository.
Any code that you would like to implement in a more general fashion should go into the framework repository.

**Important:** Have a look at the [contribution guide](CONTRIBUTING.md) which provides guidelines on how to properly work with Utopia.


### Testing
Not all test groups of the main project are available in this repository.

| Identifier        | Test description   |
| ----------------- | ------------------ |
| `model_<name>`    | The C++ model tests of model with name `<name>` |
| `models`†         | The C++ and Python tests for _all_ models |
| `models_python`†‡ | All python model tests (from `python/model_tests`) |
| `all`             | All of the above. (Go make yourself a hot beverage, when invoking this.) |

_Note:_
* Identifiers marked with `†` require all models to be built (by running
  `make all`).
* Identifiers marked with `‡` do _not_ have a corresponding `build_tests_*`
  target.
* Notice that you cannot execute tests for models of the main project in this
  repository.

#### Evaluating Test Code Coverage
Code coverage is useful information when writing and evaluating tests.
The coverage percentage of the C++ code is reported via the GitLab CI pipeline.
Check the
[`README.md` in the main repository](https://gitlab.com/utopia-project/utopia/utopia#c-code-coverage)
for information on how to compile your code with coverage report flags and how
to retrieve the coverage information.


## Dependencies

The version numbers given here are those available in the testing image, based on Ubuntu 20.04.
These version numbers are _not_ enforced.


### Additional Python Dependencies
If your model requires additional Python packages (e.g., for plotting and
testing), they are added to the [`python/requirements.txt`](python/requirements.txt) file.
All packages listed there will be automatically installed into the Python
virtual environment as part of the CMake configuration.
Note that these might require some additional dependencies to be installed via
the package manager; those are listed above.

--- 

[Utopia]: https://gitlab.com/utopia-project/utopia/utopia
