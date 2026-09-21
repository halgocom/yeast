#!/bin/bash



#container con lo stesso nome utente dell'host, user root va sostituito
cat << "EOF"
=============================================================
____    ____  _______     ___           _______.___________.
\   \  /   / |   ____|   /   \         /       |           |
 \   \/   /  |  |__     /  ^  \       |   (----`---|  |----`
  \_    _/   |   __|   /  /_\  \       \   \       |  |     
    |  |     |  |____ /  _____  \  .----)   |      |  |     
    |__|     |_______/__/     \__\ |_______/       |__|      
=============================================================
                         YEAST PIPELINE
=============================================================
EOF
while :
do
    echo $'Use devmode?\n[Y]es\n[N]o' #abilitates root access to easyly modify envs,to be eliminated this is temporary and unsafe
    read -p "" devmode
    
    case "$devmode" in
        "Y")
            user_bind="root"
        ;;
        "N")
            user_bind="$(id -u):$(id -g)"
        ;;
        *) echo "invalid option $run_choice"
        ;;
        
    esac

    read -p "Module name": module
    docker run --rm --user $user_bind --env-file ./app/modules/$module/app/env/.env -v ./app/files:/home/$module/shared_files -v ./app/modules/$module:/home/$module -it $module
done
