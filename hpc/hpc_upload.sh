#!/bin/bash
echo "Cineca Upload Utility"
read -p "Absolute or relative path of directory to upload. Type . to use the same directory as this script:": local_path;
if [[ $local_path == "." ]];then
	$local_path = $("$PWD")
fi
if [[ $local_path == ./* ]];then #checks if relative path is used
	echo $local_path
	$local_path = "$local_path" | tr "." $PWD
	echo $local_path
fi

echo $local_path
globus-url-copy -vb -cd -r $local_path   sshftp://vsunna00@gftp.g100.cineca.it:22/g100_work/IscrB_AACD/Vincenzo

