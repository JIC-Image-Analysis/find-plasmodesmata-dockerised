#!/bin/bash

CONTAINER="find-plasmodesmata-dockerised"
touch `pwd`/bash_history
docker run -it --rm -v `pwd`/bash_history:/root/.bash_history -v `pwd`/data:/data:ro -v `pwd`/scripts:/scripts:ro -v `pwd`/output:/output $CONTAINER
