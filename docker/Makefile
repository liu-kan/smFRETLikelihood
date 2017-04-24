build:
	echo "{ \n\
  "registry-mirrors": ["https://docker.mirrors.ustc.edu.cn"] \n\
  }\n\
	" > /etc/docker/daemon.json
	docker build --rm -t intelpython_mpi .
enter:
	docker exec -it mpi_node sh -c "export TERM=xterm && bash"
run:
	docker run -d --name mpi_node -p 222:22 intelpython_mpi

kill:
	-docker kill mpi_node
	-docker rm mpi_node
