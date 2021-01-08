This is a custom lightweight C++ interfact to Kodiak for use in NOVA.
This interface connects Kodiak with the NOVA using Python's ctypes library.

## Setup
1. Download [Kodiak](https://github.com/nasa/Kodiak) and all the required softwares (cmake, Boost & filib++).
2. Build a static library of Kodiak.

    Change the configuration of the Kodiak library by changing `Kodiak/CMakeLists.txt`:
    line 30 to be:
    ```
    add_definitions(-Wall -frounding-math -pedantic -Wextra -Wno-parentheses -g -O2 -fPIC)
    ```

    after that, run CMake to generate the makefiles as normal:
    ```angular2html
    cd Kodiak
    mkdir build
    cd build
    cmake ..
    ```
    then build a static library of Kodiak by:
    ```
    make KodiakStatic
    ```     
    to create the `.a` file.
    If you don't do this, the `pykodiak.so` file won't link (it will complain that -fPIC was missing in Kodiak code).
    This will produce `libkodiak.a`.
3. Now you need to create a shared library of the Kodiak interface.

    Change line 14 in `nova-public/src/utils/pykodiak/Makefile` to
    point to the directory containing `libkodiak.a`. Then run
    ```angular2html
    cd nova-public/src/utils/pykodiak
    make all
    ```
    This will produce a shared library `pykodiak.so`.

## Example
A sample usage of the Python interface can be found at `pykodiak/pykodiak_interface.py`.
    
