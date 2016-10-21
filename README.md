# find-plasmodesmata-dockerised

## Background

Scripts to identify plasmodesmata in cells.

This project is the successor of the
[find-plasmodesmata](https://github.com/JIC-CSB/find-plasmodesmata) project.


## Introduction

This image analysis project has been setup to take advantage of a technology
known as Docker.

This means that you will need to:

1. Download and install the [Docker Toolbox](https://www.docker.com/products/docker-toolbox)
2. Build a docker image

Before you can run the image analysis in a docker container.


## Build a Docker image

Before you can run your analysis you need to build your docker image.  Once you
have built the docker image you should not need to do this step again.

A docker image is basically a binary blob that contains all the dependencies
required for the analysis scripts. In other words the docker image has got no
relation to the types of images that we want to analyse, it is simply a
technology that we use to make it easier to run the analysis scripts.

```
$ cd docker
$ bash build_docker_image.sh
$ cd ..
```

## Put your images in the ``data`` directory

Put the images that you want to analyse into the ``data`` directory.


## Run the image analysis in a Docker container

The image analysis will be run in a Docker container.  The script
``run_docker_container.sh`` will drop you into an interactive Docker session.

```
$ bash run_docker_container.sh
[root@048bd4bd961c /]#
```

Now you can run the image analysis on all the images in a ``data`` directory.

```
[root@048bd4bd961c /]# python scripts/analyse_all_images.py data/ output/
```

It is also possible to run the analysis on an individual image.

```
[root@048bd4bd961c /]# python scripts/analyse_image.py data/img.lif output/
```

Or a specific series in an image.

```
[root@048bd4bd961c /]# python scripts/plasmodesmata_analysis.py data/img.lif 0 output/
```


## Specifying settings

There are three settings that can be specified. The absolute threshold used for the
segmentation (``-t`` or ``--threshold``), the minimum number of voxels required for
a spot to be classified as a plasmodesmata (``--min-voxels``) and the maximum
voxels used as an upper bound for a spot to be classified as a plasmodesmata
(``--max-voxels``).

For example to increase the maximum number of voxels allowed from 50 (the default)
to 150 one would use the command below.

```
[root@048bd4bd961c /]# python scripts/analyse_all_images.py --max-voxels 150 data/ output/
```
