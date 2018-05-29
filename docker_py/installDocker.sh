#!/bin/bash
curl -fsSL get.docker.com -o get-docker.sh
sh get-docker.sh
sudo groupadd docker
sudo usermod -aG docker $USER       
echo "!!!!! You should logout, and relogin !!!!!"