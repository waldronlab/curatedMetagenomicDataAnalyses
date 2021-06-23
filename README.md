# curatedMetagenomicAnalyses

## Docker

### Requirements

You need [Docker](https://docs.docker.com/get-docker/).

### Getting Started

First build the image:

    docker build -t "waldronlab/curatedmetagenomicanalyses" .

Then run a container based on the image:

    docker run -d -p 80:8000 -p 8787:8787 \
        --name cma waldronlab/curatedmetagenomicanalyses

Visit `0.0.0.0` in your browser to access JupyterHub. The user and password
are both `waldronlab`.

Note: RStudio is accessible at `0.0.0.0:8787`; however, it may need be
restarted to run with `docker exec -it cma rstudio-server start`.

### Start, Stop, Restart

You can start, stop, and restart the container as follows, where `cma` is the
name of the container:

    docker stop cma
    docker start cma

### Updating the image

If the image gets updated, you'll need to stop and remove the container, which
destroys any work inside the container. You must then rebuild and run the
container:

    docker stop cma
    docker rm cma
    docker build -t "waldronlab/curatedmetagenomicanalyses" .
    docker run -d -p 80:8000 -p 8787:8787 \
        --name cma waldronlab/curatedmetagenomicanalyses

If you need to copy any files from the container before you destroy it, you
can use `docker cp`:

    docker cp cma:/home/waldronlab/my_file.txt /tmp

### For testing, the following may be helpful

#### Mount a path your local system in the docker

The following mounts the current path `$(pwd)` to `/home/waldronlab/files`:

    docker run -d -p 80:8000 -p 8787:8787 \
        --volume $(pwd):/home/waldronlab/files \
        --name cma waldronlab/curatedmetagenomicanalyses

#### Execute a command in the docker container

The following will start a `bash` terminal in the container as `root`, which you
can quit using `exit`.

    docker exec -it cma bash
