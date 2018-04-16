#!/bin/bash
echo "Wait the loading finish, the console will info you ..."
docker load < smfretlikelihood.tar
docker volume create smfretlikelihood-vol
docker run -ti --rm --mount source=smfretlikelihood-vol,target=/home/liuk/data smfretlikelihood ./installApp.sh
echo "Loading finished !!!!"
sudo cp smfretlikelihood.sh /usr/local/bin/
sudo chmod a+x /usr/local/bin/smfretlikelihood.sh