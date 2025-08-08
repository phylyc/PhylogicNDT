#!/bin/bash

# method name
METHOD_NAME=phylogicndt

# dockerhub id
DOCKERHUB_ID=phylyc
VERSION=1.1

# build and push together
docker buildx build --platform linux/amd64 -t ${DOCKERHUB_ID}/${METHOD_NAME}:${VERSION} --push .
