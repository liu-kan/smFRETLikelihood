FROM fedora

MAINTAINER Liu K <liukan@126.com>

USER root
RUN dnf -y install python3-tkinter git python3-devel hdf5-devel hdf5 gcc lzo-devel bzip2-devel blosc-devel
RUN dnf -y install /usr/lib/rpm/redhat/redhat-hardened-cc1
RUN dnf -y install qt5-qtbase-gui
ENV USER liuk
RUN useradd -ms /bin/bash ${USER}
USER ${USER}
WORKDIR /home/${USER}
RUN pip3 install -i http://mirrors.aliyun.com/pypi/simple/ --trusted-host mirrors.aliyun.com --user Cython tables numpy scipy matplotlib pandas sklearn PyQt5 
RUN pip3 install -i http://mirrors.aliyun.com/pypi/simple/ --trusted-host mirrors.aliyun.com --user  lmfit notebook ipython ipywidgets
RUN git clone --depth=1 https://github.com/tritemio/FRETBursts.git
WORKDIR FRETBursts
RUN sed -i -e 's/print(cit)/#print(cit)/g' fretbursts/__init__.py
RUN python3 setup.py install --user
WORKDIR ..
RUN rm -rf FRETBursts

ADD hi.py /home/${USER}
ADD scripts/updateApp.sh /home/${USER}
ADD scripts/ptu2hf5_burstBin.sh /home/${USER}
RUN mkdir /home/${USER}/data;mkdir /home/${USER}/labdata
VOLUME /home/${USER}/data
VOLUME /home/${USER}/labdata
# USER root
# RUN dnf -y install libglvnd-glx
# RUN dnf -y install 
# USER ${USER}
# ENTRYPOINT ["/bin/bash"]
CMD ["python3", "hi.py"]

