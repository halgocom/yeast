read -p "Container name:": container
#echo "compressing folder"
#tar czvf ./CINECA/container/$container.tar.gz ./$container
#tar czvf ./CINECA/container/scripts.tar.gz ./scripts
read -p "Username:": user
read -p "Choose cluser:": cluster
echo "Sending container img"
rsync -PravzHS $PWD/singularity/CINECA/container/$container.tar.gz $user@data.$cluster.cineca.it:/g100_work/IscrB_HPCCTA/vincenzo/data/$container.tar.gz 
echo "Sending SLURM script"
rsync -PravzHS $PWD/singularity/CINECA/container/scripts.tar.gz $user@data.$cluster.cineca.it:$WORK/vincenzo/data/scripts.tar.gz 