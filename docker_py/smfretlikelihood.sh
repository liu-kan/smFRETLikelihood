#!/bin/bash

docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v $1:/home/liuk/labdata \
    --mount source=smfretlikelihood-vol,target=/home/liuk/data smfretlikelihood ./ptu2hf5_burstBin.sh $2 $3


    docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v /tmp/sm:/home/liuk/labdata \
    --mount source=smfretlikelihood-vol,target=/home/liuk/data smfretlikelihood ./ptu2hf5_burstBin.sh RSV89C224C.ptu RSV89C224C.mat