build:
	cp ~/.ssh/id_rsa.pub id_rsa.pub_
	docker build --rm -t intelpython_mpi .
enter:
	docker exec -it mpi_node bash
run:
	docker run -d --name mpi_node -v `pwd`/../:/smfret -p 22:22 -p 1000:1000 intelpython_mpi
	#docker run -d --name mpi_node -v `pwd`/../:/smfret -p 10000-10100:10000-10100 -p 222:22 intelpython_mpi

kill:
	-docker kill mpi_node
	-docker rm mpi_node

mirrordockerhub:
	echo "{ \n\
	\"registry-mirrors\": [\"https://docker.mirrors.ustc.edu.cn\"] \n\
	}\n\
	" > /etc/docker/daemon.json
