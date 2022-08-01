FROM ubuntu:focal
MAINTAINER Edward Kim EMAIL ehkim@cs.unc.edu
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt update && apt install -y cmake git g++ python3 libboost-all-dev wget python3-pip

RUN git clone https://github.com/Tarheel-Formal-Methods/kaa-dynamic
RUN git clone https://github.com/nasa/Kodiak
RUN wget http://www2.math.uni-wuppertal.de/wrswt/software/filib++/filibsrc-3.0.2.tar.gz && tar zxvf filibsrc-3.0.2.tar.gz

RUN cd filibsrc && ./configure CFLAGS=-fPIC CPPFLAGS=-fPIC CXXFLAGS=-fPIC && make install
WORKDIR Kodiak
RUN sed -i '30s/.*/add_definitions(-Wall -frounding-math -pedantic -Wextra -Wno-parentheses -g -O2 -fPIC)/' CMakeLists.txt
RUN mkdir build && cd build && cmake .. && make KodiakStatic

WORKDIR /kaa-dynamic
RUN git pull origin arch22 && git checkout -b arch22
RUN cd kaa/pykodiak && sed -i '11s|.*|KODIAK_LIB=/Kodiak/build/libkodiak.a|' Makefile && sed -i '14s|.*|KODIAK_INCLUDE=/Kodiak/src|' Makefile && make all
RUN pip3 install -r requirements.txt
