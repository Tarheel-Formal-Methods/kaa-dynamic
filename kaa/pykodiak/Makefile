# -*- indent-tabs-mode:t; -*-

# set this to the path to the .a file
# important note: you need to change the configuration of the Kodiak library by editing CMakeLists.txt
# change line 30 to be:
#
# add_definitions(-Wall -frounding-math -pedantic -Wextra -Wno-parentheses -g -O2 -fPIC)
#
# after that, you can use cmake to generate the makefiles as normal, then do 'make KodiakStatic' to create the .a file
# If you don't do this, the pykodiak.so file won't link (it will complain that -fPIC was missing in Kodiak code)
KODIAK_LIB= /Users/edwardkim/Work/Kodiak/build/libkodiak.a

# set this to the path to the kodiak header files
KODIAK_INCLUDE= /Users/edwardkim/Work/Kodiak/src

OUTPUT=pykodiak.so

CC=g++

CFLAGS = -I$(KODIAK_INCLUDE) -Wall -fPIC -frounding-math -O2 -std=c++14

LDFLAGS = -shared -lprim -lboost_serialization

all: .all

.PHONY: .all .check_kodiak_a .check_kodiak_include

.check_kodiak_a:
	@if [ ! -f "$(KODIAK_LIB)" ] ; then \
        echo "Kodiak .a file not found at path '$(KODIAK_LIB)'. Did you set KODIAK_LIB in the Makefile?"; \
        exit 1; \
	else exit 0; \
    fi

.check_kodiak_include:
	@if [ ! -f "$(KODIAK_INCLUDE)/kodiak.hpp" ] ; then \
        echo "Kodiak header file not found at path '$(KODIAK_INCLUDE)/kodiak.hpp'. Did you set KODIAK_INCLUDE in the Makefile?"; \
        exit 1; \
	else exit 0; \
    fi

$(OUTPUT): $(KODIAK_LIB) pykodiak.cpp
	$(CC) $(CFLAGS) pykodiak.cpp $(KODIAK_LIB) $(LDFLAGS) -o $(OUTPUT)

.all: .check_kodiak_a .check_kodiak_include $(OUTPUT)
	@exit 0

clean:
	rm -fv $(OUTPUT)
