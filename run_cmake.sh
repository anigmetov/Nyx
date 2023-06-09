#!/bin/bash

cmake -DCMAKE_BUILD_TYPE=Release -DNyx_TORCH=ON -DNyx_OMP=OFF -DNyx_HEATCOOL=ON -DNyx_SINGLE_PRECISION_PARTICLES=OFF -DTorch_DIR=/home/narn/code/libtorch_cuda12/share/cmake/Torch ..
