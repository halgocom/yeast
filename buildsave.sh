#!/bin/bash

./build.sh && docker image save screener:latest > ./docker/screener.tar