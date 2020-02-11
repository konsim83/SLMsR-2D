# The Semi-Lagrangian Multiscale Reconstruction Method (SLMsR)

This is the prototypical Julia implementation of the semil-Lagrangian basis reconstruction method *SLMsR*.

The aim is to demonstrate for the passive tracer transport problem with diffusion that it is
possible to reconstruct subgrid-features that can not be resolved 
in coarse-grid simulations locally into basis functions even when the
problem is dominated by advection. Standard multiscale methods for often fail in
this case.


The SLMsR code requires:

* A Linux distribution (we used Debian)
* **[Julia](www.julialang.com)** v1.3	
* **[Paraview](https://www.paraview.org/)** for the visualization

To add all necessary Julia packages just run the <code>startup.jl</code> script 
in the root folder:

```
include("startup.jl")
```


### How to obtain the code?

CaSIS requires:

* A Linux distribution (we used Debian and Ubuntu for the development)
* **cmake** v2.8.12 or higher	
* **clang-format-6.0** (available in most linux distributions)
* A working **debugger** (we use gdb)
* **doxygen**	
* A working installation of **[deal.ii](www.dealii.org)** v9.1.1 or higher 
with **MPI** and all other **Trilinos** dependencies must be installed. This
can easily be done through the **[spack](https://spack.readthedocs.io/en/latest/)** 
package manager
	
First you have to clone the project through ssh

``` 
git clone git@gitlab.rrz.uni-hamburg.de:BXXXXXX/Casis.git casis_code
```

or through https:

```
git clone https://gitlab.rrz.uni-hamburg.de/BXXXXXX/Casis.git casis_code
```
*BXXXXXX* represents you B-account number.

Note that for ssh you need a deploy key (ask a maintainer). Then we create 
the Makefile and eclipse project files using cmake. We do this in a separate 
folder (*out-of-source-build*)

```
mkdir casis_build
cd casis_build
cmake -DDEAL_II_DIR=/path/to/dealii-9.1.1/ -DCMAKE_ECLIPSE_MAKE_ARGUMENTS=-jN -G"Eclipse CDT4 - Unix Makefiles" ../casis_code
```

where N is the number of cores. This way one can import the project into
eclipse as an existing project.

```
make debug
```

switches the compile mode to debug mode (but does not compile anything) and

```
make release
```

switches to optimized mode. To compile the project on N cores you can type

```
make -jN
```

After everything is compiled you can run the code

```
mpirun -n N source/main_executable [aditional options]
```

where N is the number of MPI ranks (processes).


### How to extend the code?

The master branch is protected so each developer needs to
create a feature branch and when ready submit a merge request
to the master branch. The code can then be reviewed and fixed
if necessary. All commits should be squashed after the review
is done. An (interactive) rebase may be necessary if the master
branch has moved forward.     


### Unit Tests

The Project comes with a unit test suite based on deal.ii's test suite.
Each test is a *.cc file in the test folder under the casis_code directory.
For how to set up a test see the documentation [here](https://www.dealii.org/current/users/testsuite.html).
A specific test can be run through

```
ctest -V -R "test/name_of_test"
```

or similarly through

```
ctest -V -R
```


### Code Indentation

The code should be indented. We use the style that is also used by deal.ii
to make code review easier. The code is automatically indented through

```
make indent
```


### Documentation

A documention should be provided for all namespaces, classes, class members and
functions in the *.h file. Additionally code comments are highly appreciated. For 
the doxygen style see an example file and the [doxygen documentation](http://www.doxygen.nl/manual/index.html).
By convention doxygen comments start with 

```
/*!
 * @brief This is a brief description.
 *
 * The detailed documentation goes here.
 */ 
```

Doxygen should automatically pick up input arguments (if enable in the
eclipse project settings).