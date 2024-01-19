#!/bin/bash
build_manylinux() {
    PLAT=$1
    DOCKER_IMAGE=$2
    PRE_CMD=$3
    sudo docker run --rm -e PLAT=$PLAT -v `pwd`:/io $DOCKER_IMAGE $PRE_CMD \
        /io/tools/build-wheels-manylinux.sh
}

build_manylinux manylinux_2_28_x86_64 quay.io/pypa/manylinux_2_28_x86_64
