#!/bin/bash

TAG=tsp:latest

xhost +local:docker

has_nvidia_support() {
    docker info | grep -i nvidia > /dev/null
    return $?
}
NVIDIA_DOCKER_FLAGS=
if has_nvidia_support; then
    NVIDIA_DOCKER_FLAGS+="--gpus all                         "
    NVIDIA_DOCKER_FLAGS+="-e NVIDIA_DRIVER_CAPABILITIES=all  "
    NVIDIA_DOCKER_FLAGS+="-e __NV_PRIME_RENDER_OFFLOAD=1     "
    NVIDIA_DOCKER_FLAGS+="-e __GLX_VENDOR_LIBRARY_NAME=nvidia"
fi

# Useful resources for utilizing the Nvidia GPU card for rendering
# https://download.nvidia.com/XFree86/Linux-x86_64/535.183.01/README/primerenderoffload.html

docker run                                  \
        -it                                 \
        --rm                                \
        -u $(id -u):$(id -g)                \
        -w /ws                              \
        -v $(pwd):/ws                       \
        -e DISPLAY=$DISPLAY                 \
        -v /tmp/.X11-unix:/tmp/.X11-unix    \
        --privileged                        \
        ${NVIDIA_DOCKER_FLAGS}              \
        $TAG                                \
        "$@"
