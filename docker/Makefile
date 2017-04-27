build:
	cp ~/.ssh/id_rsa.pub id_rsa.pub_
	docker build --rm -t intelpython_mpi -f Dockerfile.mpich .
enter:
	docker exec -it mpi_node bash
run:
	docker run -d --name mpi_node -v `pwd`/../:/smfret -p 22:22 -p 1000:1100 intelpython_mpi
	#docker run -d --name mpi_node -v `pwd`/../:/smfret -p 10000-10100:10000-10100 -p 222:22 intelpython_mpi

kill:
	-docker kill mpi_node
	-docker rm mpi_node

mirrordockerhub:
	echo "{ \n\
	\"registry-mirrors\": [\"https://docker.mirrors.ustc.edu.cn\"] \n\
	}\n\
	" > /etc/docker/daemon.json

mpiexec:
	docker exec --privileged mpi_node mpirun --allow-run-as-root -n 8 -H ma,sa python3 /smfret/algo/maxLikDockerMPI.py
