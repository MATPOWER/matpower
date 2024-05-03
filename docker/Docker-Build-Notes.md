Docker Build Notes
==================

Octave Images ([`matpower/octave`][1])
--------------------------------------

When a new version _x.y.z_ of the GNU Octave Docker image
([`gnuoctave/octave`][2]) becomes available ...
- Update `VER` in `Dockerfile-octave` to default to new version.
- Add corresponding new line in [`build_octave_images.sh`][3].
- Build new version.
  ```
  docker build -f docker/Dockerfile-octave -t matpower/octave:<x.y.z> .
  ```
  _(Note: Likely requires an explicit `--platform=linux/amd64` flag when
  building on macOS with Apple Silicon chip.)_
- Update `latest` tag.
  ```
  docker tag matpower/octave:<x.y.z> matpower/octave:latest
  ```
- Push to [Docker Hub][4].
  ```
  docker push matpower/octave:<x.y.z>
  docker push matpower/octave:latest
  ```
- Rebuild `<X.Y>c` versions of MATPOWER images (including legacy), test them,
  then push them to Docker Hub.
- Update value listed for _current release_ of Octave in **Versions** section
  of [MATPOWER-Docker.md][5].
- Possibly add row to table in [Octave-Docker.md][6].
- Manually update README files for [`matpower/octave`][1] and
  [`matpower/matpower`][7] on DockerHub from [Octave-Docker.md][6] and
  [MATPOWER-Docker.md][5], respectively.


MATPOWER Images ([`matpower/matpower`][7]) (v7.0+)
------------------------------------------------

When a new MATPOWER version _X.Y_ is released ...
- Update value listed for _current release_ of MATPOWER in **Versions** section
  of [MATPOWER-Docker.md][5], and add new rows to table.
- Build the new image with labels `<X.Y>`, `<X.Y>c`, `latest` 
  ```
  docker build --build-arg BRANCH=<X.Y> -f docker/Dockerfile -t matpower/matpower:<X.Y> .
  ```
- Update `latest` tag and create `<X.Y>c` tag
  ```
  docker tag matpower/matpower:<X.Y> matpower/matpower:latest
  docker tag matpower/matpower:<X.Y> matpower/matpower:<X.Y>c
  ```

To build a new MATPOWER image from the latest `master` branch ...
```
docker build --build-arg MP_SRC=github -f docker/Dockerfile -t matpower/matpower:dev-<YYYY-MM-DD> .
```

To build for a specific MATPOWER version (_X.Y_), specify the corresponding GitHub tag
in `BRANCH` arg.
```
docker build --build-arg MP_SRC=github --build-arg BRANCH=<X.Y> -f docker/Dockerfile -t matpower/matpower:<X.Y>c .
```

To build for a specific Octave version (_x.y.z_), specify the tag in `BASE_TAG` arg.
```
docker build --build-arg MP_SRC=github --build-arg BASE_TAG=<x.y.z> --build-arg BRANCH=<X.Y> -f docker/Dockerfile -t matpower:<X.Y>-oct-<x.y.z> .
```

See [`build_matpower_images.sh`][8] for examples.

MATPOWER Legacy Images ([`matpower/matpower`][7]) (pre v7.0)
------------------------------------------------------------

See [`build_matpower_images.sh`][8] for examples.

---

[1]: https://hub.docker.com/r/matpower/octave
[2]: https://hub.docker.com/r/gnuoctave/octave
[3]: ./build_octave_images.sh
[4]: https://hub.docker.com
[5]: ./MATPOWER-Docker.md
[6]: ./Octave-Docker.md
[7]: https://hub.docker.com/r/matpower/matpower
[8]: ./build_matpower_images.sh
