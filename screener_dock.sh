#! /bin/bash

read -p "Protocol name:": protocol
micromamba run -n Screener python /home/screener/scripts/base/$protocol.py