#!/bin/bash

cat << "EOF"
==========================================
██╗   ██╗███████╗ █████╗ ███████╗████████╗
██║   ██║██╔════╝██╔══██╗██╔════╝╚══██╔══╝
██║   ██║█████╗  ███████║███████╗   ██║   
╚██╗ ██╔╝██╔══╝  ██╔══██║╚════██║   ██║   
 ╚████╔╝ ███████╗██║  ██║███████║   ██║   
  ╚═══╝  ╚══════╝╚═╝  ╚═╝╚══════╝   ╚═╝   
==========================================
            YEAST PIPELINE
==========================================
EOF

while :
do
    echo $'Want to clean build cache?\n[Y]es\n[N]o'
    read -p "" clean_choice

    case "$clean_choice" in
        "Y")
            docker system prune -a;
            break;;
        "y")
            echo "must be capitalized";;
        "N")
            break;;
        *) echo "invalid option $clean_choice,either Y or N is allowed"
        ;;
    esac
done



echo $'Insert the name of the module to build'
read -p "Module name": module
docker build -f ./app/modules/$module/Dockerfile -t $module ./app/modules/$module;

