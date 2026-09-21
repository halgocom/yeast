#! /bin/bash

read -p "Protocol name:": protocol
read -p "Options:": options
micromamba run -n crest python /home/crest/app/scripts/$protocol.py $options