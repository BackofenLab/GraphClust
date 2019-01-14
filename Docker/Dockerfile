# GraphClust Docker
#
# VERSION 1.0

# FROM ubuntu:14.04
FROM debian:wheezy

MAINTAINER Milad Miladi <miladim@informatik.uni-freiburg.de>

# =========================================================
# Setting up sge recepie and dependency packages

RUN DEBIAN_FRONTEND=noninteractive apt-get  update

RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gridengine-master gridengine-exec gridengine-common gridengine-qmon gridengine-client
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y ghostscript wget graphicsmagick graphicsmagick-imagemagick-compat r-base-core\
    cpanminus make sudo nano less vim procps apt-utils apt-transport-https bzip2 time
RUN DEBIAN_FRONTEND=noninteractive apt-get  update

RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get install --no-install-recommends -y build-essential libatlas-dev g++ 

# Usually missing in base libraries and needed by GC
RUN cpanm Math::Round 
RUN pip install -q numpy --upgrade
RUN pip install -q seaborn
RUN pip install -q six

RUN mkdir -p /usr/local/user/

RUN useradd bisge001 --create-home
RUN adduser bisge001 sudo

ADD exec_host /exec_host
# Set maximum of available memory - this is avaiable memory - 1GB
RUN sed -i "s|complex_values        NONE|complex_values        h_vmem=`grep 'MemTotal' /proc/meminfo | awk '{print ($2 - 1000000)}'`k|g" /exec_host

ADD host_group_entry /host_group_entry
ADD queue /queue

USER root
# Expose port 8080 to the host
EXPOSE :8080


# =========================================================
# setting up GraphClust

# Set necessary ENV variables
ENV HOST_URL="ftp://biftp.informatik.uni-freiburg.de/pub/RNAscClust/"\
    PATH="/usr/local/perl/5.12.2/bin/":$PATH

# Below Needed for precompiled PETfold
RUN echo 'deb http://ftp.debian.org/debian experimental main' >> /etc/apt/sources.list &&\
    apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y -t experimental libc6-dev


# GraphClust and dependency tools (built binaries)
# tools-binaries created with this command on bi ui server: tar -czf binaries.tar.gz "/usr/local/user/GraphClust-0.7.6/" "/usr/local/user/locarna/1.7.16/" "/usr/local/rnashapes/2.1.6/" "/usr/local/user/ViennaRNA/2.1.7/" "/usr/local/user/ViennaRNA/2.1.2/" "/usr/local/infernal/1.0.2/" "/usr/local/rnaz/2.1-a0130d9/" "/usr/local/ncbiblast/2.2.15/" "/usr/local/user/cmfinder-0.2/" "/usr/local/user/PETfold/"
# 
# perl 5.12 (pre built) Todo: install perl
RUN wget -qO-  "$HOST_URL/tools-binaries.tar.gz" | tar -zxf - -C / &&  \
    wget -qO- "$HOST_URL/perl-512.tar.gz" | tar -zxf - -C / 
   

# resolve tty error msg in running GraphClust sge run as root
RUN sed -i -- 's/mesg n/tty -s \&\& mesg n/g' /root/.profile

# RNAscClust source and configs

        

# startup commands to setup the server and etc
ADD run.sh /root/run.sh
RUN chmod +x /root/run.sh

CMD bash -C '/root/run.sh';'bash'
