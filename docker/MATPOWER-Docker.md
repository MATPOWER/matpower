MATPOWER Docker Image
=====================

A complete MATPOWER desktop environment is available in a [Docker][1]
container by using the [matpower/matpower-desktop][2] image on [Docker Hub][3].
This image essentially just adds [MATPOWER][4] and the [MATPOWER Extras][5]
to the [compdatasci/octave-desktop][6] image, which provides an [Ubuntu][7]
desktop environment with X Windows and [GNU Octave][8], displaying in your
web-browser.

It can also be used in command-line mode without the graphical user
interface (GUI), if that is preferred.


System Requirements
-------------------

You will need working installations of:
- [Docker][9], and
- [Python][10] _(optional, required for use of the GUI)_


Using the MATPOWER Docker Image
-------------------------------

To set up and use the MATPOWER Docker image, please follow the instructions
found in the [README][11] for the [compdatasci/octave-desktop][6] Docker image,
with the following minor difference:

- The Python script is called [`matpower_desktop.py`][12] instead of
  `octave_desktop.py` and is found at
  [`https://raw.githubusercontent.com/MATPOWER/matpower/master/docker/matpower_desktop.py`][12]
  (not at `https://raw.githubusercontent.com/compdatasci/octave-desktop/master/octave_desktop.py`).


Quick Start
-----------

#### Run via Python script — Access GUI through Web Browser
With Python installed, Docker installed and running, and the
[`matpower_desktop.py`][12] script in your working directory ...
```
python matpower_desktop.py -p
```

Then via the Terminal window of the Ubuntu desktop that appears in your
web-browser window ...
```
octave &
```

This runs Octave in your newly launched Docker container, with [MATPOWER][4]
and the [MATPOWER Extras][5] already installed in the Octave path.

You can also access Octave and [MATPOWER][4] in this container (the one
you launched via Python script) from the command-line on your host machine,
either via SSH ...
```
ssh -X -p 2222 ubuntu@localhost
```
... then, in the container ...
```
octave-cli
```

... or directly via `docker exec` ...
```
docker exec -it <container-name> octave-cli
```
... where `<container-name>` can be found via `docker container ls --all`.


#### Run Directly — Command-line Only

If you do not have a need for the GUI, you can launch a shell or command-line
Octave session directly in a new container (without the Python script) via
`docker run`. The following command runs a new container from the MATPOWER
Docker image, leaving you at the Octave command-line.
```
docker run --rm -it matpower/matpower-desktop:latest /sbin/my_init -- octave-cli
```
Or, if you prefer to enter a session where you can use shell commands
(including `octave-cli`), use:
```
docker run --rm -it matpower/matpower-desktop:latest /sbin/my_init -- bash -l
```
In either case, [MATPOWER][4] and the [MATPOWER Extras][5] are already
installed in the Octave path.

---

[1]: https://www.docker.com
[2]: https://hub.docker.com/r/matpower/matpower-desktop
[3]: https://hub.docker.com/
[4]: https://github.com/MATPOWER/matpower
[5]: https://github.com/MATPOWER/matpower-extras
[6]: https://hub.docker.com/r/compdatasci/octave-desktop
[7]: https://www.ubuntu.com
[8]: https://octave.org
[9]: https://www.docker.com/products/docker-desktop
[10]: https://www.python.org
[11]: https://github.com/compdatasci/octave-desktop/blob/master/README.md
[12]: https://raw.githubusercontent.com/MATPOWER/matpower/master/docker/matpower_desktop.py
