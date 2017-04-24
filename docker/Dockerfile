FROM intelpython/intelpython3_full

MAINTAINER Liu K <liukan@126.com>

USER root

RUN echo -e "deb http://mirrors.aliyun.com/debian jessie main contrib non-free\n\
deb-src http://mirrors.aliyun.com/debian jessie main contrib non-free\n\
deb http://mirrors.aliyun.com/debian jessie-proposed-updates main contrib non-free\n\
deb-src http://mirrors.aliyun.com/debian jessie-proposed-updates main contrib non-free\n\
deb http://mirrors.aliyun.com/debian jessie-updates main contrib non-free\n\
deb-src http://mirrors.aliyun.com/debian jessie-updates main contrib non-free\n\
deb http://mirrors.ustc.edu.cn/debian-security jessie/updates main\n\
" > /etc/apt/sources.list

RUN apt-get -y update && \
    apt-get -y vim mpich openssh-server && \
    conda install -y pyqt

## sshd
RUN mkdir /var/run/sshd
RUN echo 'root:screencast' | chpasswd
RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile
EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
