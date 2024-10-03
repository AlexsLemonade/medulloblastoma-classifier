# script to get docker image up and running

sudo snap install docker

sleep 3 && sudo chmod 666 /var/run/docker.sock

docker pull ghcr.io/alexslemonade/medulloblastoma-classifier:latest

docker run -d --rm --mount type=volume,dst=/home/rstudio,volume-driver=local,volume-opt=type=none,volume-opt=o=bind,volume-opt=device=/home/ubuntu -e PASSWORD=onecupatatime ghcr.io/alexslemonade/medulloblastoma-classifier

container_id=$(docker ps -ql)

docker exec -it $container_id bash

