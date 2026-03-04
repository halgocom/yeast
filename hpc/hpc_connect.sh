#!/bin/bash
echo "Generating certificates..."
step ca bootstrap --ca-url=https://sshproxy.hpc.cineca.it --fingerprint 2ae1543202304d3f434bdc1a2c92eff2cd2b02110206ef06317e70c1c1735ecd
eval $(ssh-agent)
read -p "email:": email;
step ssh login $email  --provisioner cineca-hpc
read -p "Username:": user
read -p "Choose cluser:": cluster
echo "Connecting to cluster..."
ssh-keygen -f "$HOME/.ssh/known_hosts" -R "login.$cluster.cineca.it"
ssh $user@login.$cluster.cineca.it -o hashknownhosts=no

