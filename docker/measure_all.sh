#/bin/sh

rm -Rf results
mkdir results
sudo docker build -t kaa-arch .
sudo docker run -v $PWD/results:/kaa-dynamic/figures -v $PWD/results:/kaa-dynamic/data kaa-arch:latest \
      python3 -c "from examples.test_arch22 import *; run_all_benchmarks()"
