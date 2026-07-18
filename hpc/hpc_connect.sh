#!/bin/bash
echo "Generating certificates..."
step ca bootstrap --ca-url=https://sshproxy.hpc.cineca.it --fingerprint MY_PRINT
eval $(ssh-agent)
read -p "email:": email;
step ssh login $email  --provisioner cineca-hpc
read -p "Username:": user
read -p "Choose cluser:": cluster
echo "Connecting to cluster..."
ssh-keygen -f "$HOME/.ssh/known_hosts" -R "login.$cluster.cineca.it"
ssh $user@login.$cluster.cineca.it -o hashknownhosts=no

