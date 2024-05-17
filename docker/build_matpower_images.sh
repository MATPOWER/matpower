#!/bin/sh

## Dockerfile-legacy
## MATPOWER images for versions 4.0, 4.1, 5.0, 5.1, 6.0
## with ~ contemporary Octave version
docker build --build-arg VER=4.0 --build-arg BASE_TAG=4.0.3 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:4.0 .
docker build --build-arg VER=4.1 --build-arg BASE_TAG=4.0.3 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:4.1 .
docker build --build-arg VER=5.0 --build-arg BASE_TAG=4.0.3 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.0 .
docker build --build-arg VER=5.1 --build-arg BASE_TAG=4.0.3 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.1 .
docker build --build-arg VER=6.0 --build-arg BASE_TAG=4.2.2 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:6.0 .

## with current (latest) Octave version
docker build --build-arg VER=4.0 --build-arg WARN=1 -f docker/Dockerfile-legacy -t matpower/matpower:4.0c .
docker build --build-arg VER=4.1 --build-arg WARN=1 -f docker/Dockerfile-legacy -t matpower/matpower:4.1c .
docker build --build-arg VER=5.0 --build-arg WARN=1 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.0c .
docker build --build-arg VER=5.1 --build-arg WARN=1 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.1c .
docker build --build-arg VER=6.0 --build-arg WARN=1 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:6.0c .

# ## other specific Octave versions
# docker build --build-arg VER=4.0 --build-arg BASE_TAG=4.2.2 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:4.0-oct-4.2.2 .
# docker build --build-arg VER=4.0 --build-arg BASE_TAG=4.4.1 -f docker/Dockerfile-legacy -t matpower/matpower:4.0-oct-4.4.1 .
# docker build --build-arg VER=4.0 --build-arg BASE_TAG=5.2.0 -f docker/Dockerfile-legacy -t matpower/matpower:4.0-oct-5.2.0 .
# docker build --build-arg VER=4.1 --build-arg BASE_TAG=4.2.2 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:4.1-oct-4.2.2 .
# docker build --build-arg VER=4.1 --build-arg BASE_TAG=4.4.1 -f docker/Dockerfile-legacy -t matpower/matpower:4.1-oct-4.4.1 .
# docker build --build-arg VER=4.1 --build-arg BASE_TAG=5.2.0 -f docker/Dockerfile-legacy -t matpower/matpower:4.1-oct-5.2.0 .
# docker build --build-arg VER=5.1 --build-arg BASE_TAG=4.2.2 --build-arg MORE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.1-oct-4.2.2 .
# docker build --build-arg VER=5.1 --build-arg BASE_TAG=4.4.1 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.1-oct-4.4.1 .
# docker build --build-arg VER=5.1 --build-arg BASE_TAG=5.2.0 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:5.1-oct-5.2.0 .
# docker build --build-arg VER=6.0 --build-arg BASE_TAG=4.4.1 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:6.0-oct-4.4.1 .
# docker build --build-arg VER=6.0 --build-arg BASE_TAG=5.2.0 --build-arg PATCH_PSSE=1 -f docker/Dockerfile-legacy -t matpower/matpower:6.0-oct-5.2.0 .

## Dockerfile
## MATPOWER images for version 7.0 and newer
## with contemporary Octave version
docker build --build-arg MP_SRC=github --build-arg BRANCH=7.0 --build-arg BASE_TAG=5.1.0 -f docker/Dockerfile -t matpower/matpower:7.0 .
docker build --build-arg MP_SRC=github --build-arg BRANCH=7.1 --build-arg BASE_TAG=5.2.0 -f docker/Dockerfile -t matpower/matpower:7.1 .
docker build --build-arg MP_SRC=github --build-arg BRANCH=8.0b1 --build-arg BASE_TAG=7.3.0 -f docker/Dockerfile -t matpower/matpower:8.0b1 .
docker build --build-arg MP_SRC=github --build-arg BRANCH=8.0 --build-arg BASE_TAG=9.1.0 -f docker/Dockerfile -t matpower/matpower:8.0 .
## with current (latest) Octave version
docker build --build-arg MP_SRC=github --build-arg BRANCH=7.0 --build-arg WARN=1 -f docker/Dockerfile -t matpower/matpower:7.0c .
docker build --build-arg MP_SRC=github --build-arg BRANCH=7.1 --build-arg PATCH_MOST_TEST=1 -f docker/Dockerfile -t matpower/matpower:7.1c .
docker build --build-arg MP_SRC=github --build-arg BRANCH=8.0b1 --build-arg PATCH_MOST_TEST=1 -f docker/Dockerfile -t matpower/matpower:8.0b1c .
docker build --build-arg MP_SRC=github --build-arg BRANCH=8.0 --build-arg PATCH_MOST_TEST=1 -f docker/Dockerfile -t matpower/matpower:8.0c .
# docker build --build-arg MP_SRC=github --build-arg BRANCH=7.1 --build-arg PATCH_MOST_TEST=1 -f docker/Dockerfile -t matpower/matpower:latest .
docker tag matpower/matpower:8.0c matpower/matpower:latest

## dev versions from current master branch
docker build --build-arg MP_SRC=github --build-arg BASE_TAG=4.4.1 -f docker/Dockerfile -t matpower/matpower:dev-latest-4.4.1 .
docker build --build-arg MP_SRC=github --build-arg BASE_TAG=5.2.0 -f docker/Dockerfile -t matpower/matpower:dev-latest-5.2.0 .
docker build --build-arg MP_SRC=github -f docker/Dockerfile -t matpower/matpower:dev-2022-01-26 .
docker build --build-arg MP_SRC=github -f docker/Dockerfile -t matpower/matpower:dev-latest .

## local dev build
# docker build -f docker/Dockerfile -t matpower/matpower:dev-local .

# MATPOWER version - latest Octave version at release - min required Octave ver
# 4.0 - 3.2.4 - 3.2+    use matpower/octave:4.0.3
# 4.1 - 3.4.3 - 3.2+    use matpower/octave:4.0.3
# 5.0 - 3.8.2 - 3.4+    use matpower/octave:4.0.3
# 5.1 - 3.8.2 - 3.4+    use matpower/octave:4.0.3
# 6.0 - 4.2.0 - 3.4+    use matpower/octave:4.2.2
# 7.0 - 5.1.0 - 4.0+    use matpower/octave:5.1.0
# 7.1 - 5.2.0 - 4.0+    use matpower/octave:5.2.0
# 
# # ## required by MATPOWER 5.0 with Octave x.x
# warning('off', 'Octave:num-to-str');
# warning('off', 'Octave:data-file-in-path');
# 
# ## required by MATPOWER 5.x/6.x with Octave >= 4.4+
# cd $HOME/packages/matpower
# patch -u -b psse_parse_line.m -i psse_parse_line5.patch
# patch -u -b psse_parse_section.m -i psse_parse_section5.patch
# 
# ## required by MATPOWER 6.0 with all versions of Octave
# cd $HOME/packages/matpower
# patch -u -b have_fcn.m -i have_fcn6.patch
# 
# ## required by MATPOWER 6.0, 7.0 with Octave 6.x
# warning('off', 'Octave:empty-index');

# for debugging, may want to try:
#   DOCKER_BUILDKIT=0 docker build ...
