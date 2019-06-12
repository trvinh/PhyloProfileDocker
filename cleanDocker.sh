#!/bin/bash

# remove exited containers:
docker rm $(docker ps -qa --no-trunc --filter "status=exited")

# remove unused images:
docker rmi $(docker images --filter "dangling=true" -q --no-trunc) -f

# remove unused volumes:
# find '/var/lib/docker/volumes/' -mindepth 1 -maxdepth 1 -type d | grep -vFf <(
#   docker ps -aq | xargs docker inspect | jq -r '.[] | .Mounts | .[] | .Name | select(.)'
# ) | xargs -I rm -fr

### build image
# sudo docker build --build-arg WHEN=2019-06-03 -t phyloprofile .

### run image
# docker run -p 8080:80 phyloprofile

### tag image
# docker tag c7b3a6ffcf11 trvinh/phyloprofile:latest

### push to docker hub
# docker push trvinh/phyloprofile
