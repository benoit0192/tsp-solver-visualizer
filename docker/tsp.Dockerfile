FROM nvidia/cuda:12.2.2-devel-ubuntu22.04

RUN apt-get update                  && \
    apt-get install -y                 \
        build-essential                \
        git                            \
        cmake                          \
        valgrind                       \
        gdb                            \
        libboost-all-dev               \
        libx11-dev                     \
        libglfw3-dev                   \
        libglew-dev                    \
        mesa-common-dev                \
        mesa-utils                     \
        libglm-dev                  && \
    apt-get clean                   && \
    rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
