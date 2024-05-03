MATPOWER's Octave Docker Image
==============================

A complete [GNU Octave][1] environment is available in a [Docker][2]
container by using the [matpower/octave][3] image on [Docker Hub][4].
This image is used as base for the [MATPOWER][5] Docker image
([matpower/matpower][6]) and consists of the corresponding official [GNU
Octave][1] image ([gnuoctave/octave][7]) with MEX interfaces for several
optimization packages pre-installed.

Specifically, the following packages are included, depending on the
version of the Octave image being used.

| Octave Version | [IPOPT][8] | [OSQP][9] | [SeDuMi][10] | [SDP3][11] | [YALMIP][12] |
| :------------: | :--------: | :-------: | :----------: | :--------: | :----------: |
|     9.x.x      |  &check;   |  &check;  |   &check;    |  &check;   |   &check;    |
|     8.x.x      |  &check;   |  &check;  |   &check;    |  &check;   |   &check;    |
|     7.x.x      |  &check;   |  &check;  |   &check;    |  &check;   |   &check;    |
|     6.x.x      |  &check;   |  &check;  |   &check;    |  &check;   |   &check;    |
|     5.x.x      |  &check;   |  &check;  |   &check;    |  &check;   |   &check;    |
|     4.4.x      |  &check;   |  &check;  |   &check;    |  &check;   |   &check;    |
|     4.2.x      |  &check;   |           |   &check;    |  &check;   |   &check;    |
|     4.0.x      |            |           |   &check;    |  &check;   |   &check;    |

All packages are built from the latest versions of the source from GitHub,
except YALMIP, which uses R20180817.

This image supports both a simple command-line mode and, with an X11
server running on the host, the full Octave GUI.


System Requirements
-------------------

You will need working installations of:
- [Docker][13], and
- an X11 server _(optional, required for use of the GUI)_


Getting Started
---------------

### 1. Get the Octave Docker image
```
docker pull docker.io/matpower/octave:latest
```

### 2. Start Octave in the container

**_command line only_**
```
docker run -it --rm matpower/octave:latest octave-cli
```

**_graphical user interface_**  
_(requires X11 server to be running and `DISPLAY` environment variable
to be set properly)_
```
docker run -it --rm --network=host --env="DISPLAY" \
  --volume="$HOME/.Xauthority:/root/.Xauthority:rw" \
  matpower/octave:latest octave --force-gui
```

This runs [Octave][1] in your newly launched Docker container, with
the included packages already installed in the Octave path.


### 3. Exectute Octave commands

For example ...
```
ver
```


Additional Notes
----------------

- To run a different version of [Octave][1], replace the `latest` tag
  with the appropriate available version tag.

- You can also replace `octave-cli` or `octave --force-gui` in the
  `docker run` command with `bash` to start your container at the shell
  prompt. From there you can, for example, start multiple GUI instances of
  [Octave][1] in the same container. Two useful aliases defined in the shell
  are `ot` and `otg` to start command-line and GUI versions of [Octave][1],
  respectively.

- You can also access [Octave][1] in your running container from the
  command-line on your host machine via `docker exec` where the
  `<container-name>` can be found via `docker container ls --all`.
  ```
  docker exec -it <container-name> octave-cli
  ```


---

[1]: https://octave.org
[2]: https://www.docker.com
[3]: https://hub.docker.com/r/matpower/octave
[4]: https://hub.docker.com/
[5]: https://matpower.org
[6]: https://hub.docker.com/r/matpower/matpower
[7]: https://hub.docker.com/r/gnuoctave/octave
[8]: https://coin-or.github.io/Ipopt/
[9]: https://osqp.org
[10]: https://github.com/sqlp/sedumi
[11]: https://github.com/sqlp/sdpt3
[12]: https://yalmip.github.io
[13]: https://www.docker.com/products/docker-desktop
