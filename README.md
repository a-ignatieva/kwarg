# KwARG

Software implementing a parsimony-based heuristic algorithm for reconstructing ancestral recombination graphs (ARGs) with recurrent mutation. 

## Citation
Ignatieva, A., Lyngs&oslash;, R. B., Jenkins, P. A., and Hein, J. (2020). KwARG: Parsimonious reconstruction of ancestral recombination graphs with recurrent       mutation. [[arXiv]](http://arxiv.org/abs/2012.09562)

## License information
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Installation

This repository contains the program KwARG, and two small programs for data manipulation called "simplify" and "flip". 

A copy of Beagle is also provided, which computes exactly the minimum number of crossover recombinations required for a dataset under the infinite sites assumption (Lyngsø, R., Song, Y.S., and Hein, J. (2005). Minimum Recombination Histories by Branch and Bound. Proceedings of WABI 2005, Lecture Notes in Computer Science, 3692, pp. 239-250).

## Binaries
Binaries are available for Linux, Mac OS and Windows. 

Download the [latest release](https://github.com/a-ignatieva/kwarg/releases) for your OS and unzip. The resulting folder will contain the binaries for KwARG, simplify, flip and Beagle, as well as example input files and a tutorial explaining usage.

The programs are run via command line: through Terminal on Mac OS or Linux, or via Command Prompt on Windows. Navigate to the appropriate folder, using for example
```sh
cd Downloads/binaries_mac
```

Try the following to check if everything works correctly:
```sh
./kwarg < kreitman_snp.txt
```
on Linux and Mac OS, or
```sh
kwarg < kreitman_snp.txt
```
on Windows.

Option descriptions can be printed by running "./kwarg -h" or "kwarg -h".


## Compile from source

Clone the repository using git via command line:
```sh
git clone https://github.com/a-ignatieva/kwarg.git
```
or download the source code for the [latest release](https://github.com/a-ignatieva/kwarg/releases) and unzip.

Then follow the OS-specific instructions below to compile. Examples and usage instructions are provided in the "binaries_X" folders.

### Linux and Mac OS

Open Terminal and navigate to the "source" folder, for instance:
```sh
cd Downloads/kwarg-1.0/source
```

The following command will compile KwARG, Beagle, simplify and flip:
```sh
make all
```

The following command will delete the compiled objects:
```sh
make clean
```

### Windows

Rename "Makefile_windows" to "Makefile", and "generatebitfunctions_windows.py" to "generatebitfunctions.py" (overwriting the existing files).

Install a gcc compiler if you do not have this; for example, download WinLibs (http://winlibs.com/) to get MinGW-w64, and follow the installation instructions to set this up.
Add the mingw64/bin directory to PATH. 

You will also need python; the Makefile presumes you invoke python with the "py" command. 

Open Command Prompt from the Start menu, and navigate to the "source" folder, for instance:
```sh
cd Downloads/kwarg-1.0/source
```

Compile using 
```sh
mingw32-make all
```

# Bug reports and requests
Please leave a report using the issue tracking system, or email anastasia.ignatieva@warwick.ac.uk

# Version history

v1.0 (17/12/2020)

First release
