#!/bin/bash



#container con lo stesso nome utente dell'host, user root va sostituito 
docker run --user "$(id -u):$(id -g)" --env-file ./app/env/.env -v ./app:/home/screener -it screener
#esegue in container costruito ,rm rimuove l'istanza creata da run, v usa il volume ovvero lo spazio sul disco locale a cui il container può accedere, e quale immagine usare ovvero screener
