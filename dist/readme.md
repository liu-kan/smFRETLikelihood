# smfretlikelihood 说明文件

## 安装
解压安装包文件 smfretlikelihood-*.tar.xz
    
    tar -Jxf smfretlikelihood-0.1.tar.xz

依次执行 installDocker.sh 注销当前用户，重新登录 installApp.sh

    ./installDocker.sh 
    #logout and login, cd to smfretlikelihood path
    ./installApp.sh

## 运行

- 建立工作目录

    cd 
    mkdir smfretdata

- 将待处理的ptu文件复制到此目录

    cp LS33_RSV21c224c_alex488cy5_32MHZ1.ptu smfretdata

- 处理数据

    smfretlikelihood /home/yourname/smfretdata LS33_RSV21c224c_alex488cy5_32MHZ1.ptu fret.mat
    # /home/yourname/smfretdata 代表绝对路径

- fret.mat就是下一步可交给matlab处理的数据