# curatedMetagenomicAnalyses

## Docker

### Requirements

You need [Docker](https://docs.docker.com/get-docker/).

### Getting Started

First build the image:

    docker build -t "waldronlab/curatedmetagenomicanalyses" .

Then run a container based on the image with your password:

    docker run -d -p 80:8888 -e PASSWORD=mypass --name cma \
      waldronlab/curatedmetagenomicanalyses

Visit `localhost` in your browser. The default user and password are `waldronlab`.
