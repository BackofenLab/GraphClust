#!/usr/bin/env bash
/etc/init.d/gridengine-master start
/etc/init.d/gridengine-exec start
sudo -u sgeadmin qconf -am bisge001 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -au bisge001 users 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -Ae /exec_host 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -Ahgrp /host_group_entry 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -aattr hostgroup hostlist dockersgeserver @allhosts 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -Aq /queue 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -aattr queue hostlist @allhosts main.q 1>>qconf.log 2>>qconf.log
sudo -u bisge001 qconf -as dockersgeserver 1>>qconf.log 2>>qconf.log

# set maximum of avaiable CPU's
sed -i "s|processors            4|processors            `num=$(grep ^processor /proc/cpuinfo | wc -l) && echo $((num-1))`|g" /queue
sed -i "s|slots                 4|slots                 `num=$(grep ^processor /proc/cpuinfo | wc -l) && echo $((num-1))`|g" /queue


# CMD ["$HOME"]
mkdir '/tmp/output/'
mkdir '/RESULTS_GC/'
echo
echo "################################################################################"
echo "######################### GraphClust Docker container ##########################"
echo "################################################################################"
echo "### Run GraphClust pipeline with the sample data set cliques-low by invoking:"
echo '>>> /usr/local/user/GraphClust-0.7.6/bin/MASTER_GraphClust.pl --help'
echo 
echo "################################################################################"
echo


