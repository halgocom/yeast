#!/bin/bash

docker system prune -a
#docker build -f ./docker/Dockerfile --progress=plain -t screener:1.0 .
docker build -f ./docker/Dockerfile -t screener .
#build costruisce l'immagine e lo rende disponibile all'esecuzione
#t name of the container nome:tag
