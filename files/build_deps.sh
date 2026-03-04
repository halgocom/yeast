#! /bin/bash


#install meeko
#git clone https://github.com/forlilab/Meeko.git /home/screener/install/files/Meeko || exit
#cd /home/screener/install/files/Meeko
#git checkout develop
#micromamba run -n Screener pip install --use-pep517 /home/screener/install/files/Meeko || exit

#install scrubber
#git clone https://github.com/forlilab/scrubber.git /home/screener/install/files/scrubber || exit
#micromamba run -n Screener pip install --use-pep517 /home/screener/install/files/scrubber || exit

#compile autogrid
#git clone https://github.com/ccsb-scripps/AutoGrid /home/screener/install/files/Autogrid || exit
#cd /home/screener/install/files/Autogrid || exit
#meson setup builddir || exit
#cd /home/screener/install/files/Autogrid/builddir || exit
#meson compile || exit
#mv /home/screener/install/files/Autogrid/builddir/autogrid4 /home/screener/bin/autogrid4 || exit
#cd /home/screener/install/files

#install mordred for descriptor calculation
#micromamba run -n Screener pip install mordred || exit

#install protenix
#micromamba run -n Screener pip install protenix || exit
#git clone https://github.com/bytedance/Protenix.git /home/screener/install/files/Protenix || exit
#mv -f /home/screener/install/protenix_setup.py /home/screener/install/files/Protenix/setup.py || exit
#sed '21 i  requirements_file="/home/screener/install/files/Protenix/requirements.txt" ' /home/screener/install/files/Protenix/setup.py
#sed -i '23s/.*/with open(requirements_file) as f:/' /home/screener/install/files/Protenix/setup.py || exit
#cd /home/screener/install/files/Protenix
#micromamba run -n Screener python /home/screener/install/files/Protenix/setup.py develop  --cpu || exit

#download monomer library for cctbx
#wget -P /home/screener/install/files/monomer --no-check-certificate --content-disposition https://github.com/cctbx/cctbx_project/releases/download/v2025.10/chem_data-2025.10-pyhe0d8492_0.conda

#install monomer library for cctbx
#micromamba install -y -n Screener /home/screener/install/files/monomer/chem_data-2025.10-pyhe0d8492_0.conda || exit

#install protenix-dock
#git clone https://github.com/bytedance/Protenix-Dock.git /home/screener/install/files/Protenix_dock || exit
#mv -f /home/screener/install/pxdock_setup.py /home/screener/install/files/Protenix_dock/setup.py || exit
#sed -i '73s/.*/with open("\\home\\screener\\install\\files\\Protenix_dock\\requirements.txt") as f:/' /home/screener/install/files/Protenix_dock/setup.py || exit
#sed -i '1d' /home/screener/install/files/Protenix_dock/environment.yml || exit
#micromamba create -n pxdock -f /home/screener/install/files/Protenix_dock/environment.yml || exit
#micromamba run -n pxdock python3 /home/screener/install/files/Protenix_dock/setup.py install  || exit

#flush install files
#rm -rf /home/screener/install/files