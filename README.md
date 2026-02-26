[![](https://img.shields.io/docker/v/trvinh/phyloprofile?label=Docker%20Hub)](https://hub.docker.com/r/trvinh/phyloprofile)
[![](https://img.shields.io/docker/pulls/trvinh/phyloprofile)](https://hub.docker.com/r/trvinh/phyloprofile)

# PhyloProfileDocker
Docker image for PhyloProfile tool

# HOW-TO-USE

## Use default taxonomy database
Run the following command to download (if necessary) and start the PhyloProfile Docker container:
```
docker run --platform linux/amd64 -p 8080:80 trvinh/phyloprofile
```
Then open web browser and start PhyloProfile tool using this address:
`localhost:8080`

## Use a User-Defined Taxonomy Database

### Mount a Folder to the Container

You can mount your local taxonomy data folder into the container so that the app can use it:

```
docker run --platform linux/amd64 \
  -v /your/data/folder:/srv/shiny-server/data \
  -p 8080:80 \
  trvinh/phyloprofile
```

**NOTE:** Ensure that your mounted folder contains a `preProcessedTaxonomy.txt` file. You can generate this file using the `generated_taxonomy.R` script located in the `app` folder
```
Rscript app/generate_taxonomy.R /your/data/folder
```

### Copy Files into a Running Container

1. Find the Container ID

```
docker ps
```

Example output:
```
CONTAINER ID   IMAGE                    COMMAND                  PORTS               NAMES
b1150680c119   trvinh/phyloprofile      "/usr/bin/shiny-serv…"   0.0.0.0:8080->80    phyloprofile_app
```

2. Copy Files from Host to Container

```
docker cp /your/data/folder/. <container_id>:/srv/shiny-server/data/
```

**IMPORTANT:** Make sure `/your/data/folder` contains all (and only) the following files:
- idList.txt
- preCalcTree.nw
- rankList.txt
- taxonomyMatrix.txt
- newTaxa.txt
- taxonNamesReduced.txt

# HOW TO MAINTAIN THE DOCKER IMAGE

1. Clone or fork this repo
```
git clone https://github.com/trvinh/PhyloProfileDocker.git
```
2. Update the **app** folder reflect the latest version of (PhyloProfile)(https://github.com/BIONF/PhyloProfile)
3. Update Docker-specific files: edit the `Dockerfile` as needed (for example, to update the version of the *bioconductor/bioconductor_docker* base image)
4. (OPTIONAL) Pre-generate `preProcessedTaxonomy.txt` to save time during image build
```
Rscript app/generate_taxonomy.R app/data
```
5. Build the docker image
```
docker build . -t trvinh/phyloprofile[:1.20.4] --platform linux/amd64 --progress=plain
```

_**Note**: you need to replace <kbd>trvinh</kbd> by your Docker Username. <kbd>:0.1.2</kbd> specifies the TAG of your build (if empty, the default tag <kbd>latest</kbd> will be applied)_

6. Push to Docker Hub
```
docker push trvinh/phyloprofile[:1.20.4]
```

7. Clean build caches after pushing

Use this command to see how large is the cache
```
docker system df
```

and run this to free the disk space

```
docker buildx prune
```

__Check this [document](https://docs.docker.com/docker-hub/repos/) for more info!__

8. (OPTIONAL) Inspect the shiny server

```
docker run --platform linux/amd64 -it --rm --entrypoint /bin/bash trvinh/phyloprofile
```
