Docker Build Notes
==================

Octave Images ([`matpower/octave`][1])
--------------------------------------

When a new version x.y.z of the GNU Octave Docker image
([`gnuoctave/octave`][2]) becomes available ...
- Update `VER` in `Dockerfile-octave` to default to new version.
- Add corresponding new line in [`build_octave_images.sh`][3].
- Build new version.
  ```
  docker build -f docker/Dockerfile-octave -t matpower/octave:x.y.z .
  ```
- Update `latest` tag.
  ```
  docker tag matpower/octave:x.y.z matpower/octave:latest
  ```
- Push to [Docker Hub][4].
  ```
  docker push matpower/octave:x.y.z
  docker push matpower/octave:latest
  ```
- Rebuild `X.Yc` versions of MATPOWER images (including legacy), test them,
  then push them to Docker Hub.
- Update value listed for _current release_ of Octave in **Versions** section
  of [MATPOWER-Docker.md][5].
- Possibly add row to table in [Octave-Docker.md][6].
- Update README files for [`matpower/octave`][1] and [`matpower/matpower`][7]
  from [Octave-Docker.md][6] and [MATPOWER-Docker.md][5], respectively.


MATPOWER Images ([`matpower/matpower`][7]) (v7.0+)
------------------------------------------------

When a new MATPOWER version X.Y is released ...
- Update value listed for _current release_ of MATPOWER in **Versions** section
  of [MATPOWER-Docker.md][5], and add new rows to table.
- Build the new image with labels `X.Y`, `X.Yc`, `latest` 
  ```
  docker build --build-arg BRANCH=X.Y -f docker/Dockerfile -t matpower/matpower:X.Y .
  ```
- Update `latest` tag and create `X.Yc` tag
  ```
  docker tag matpower/matpower:X.Y matpower/matpower:latest
  docker tag matpower/matpower:X.Y matpower/matpower:X.Yc
  ```

To build a new MATPOWER image from the latest `master` branch without
[MP-Element][8] ...
```
docker build -f docker/Dockerfile -t matpower/matpower:dev-YYYY-MM-DD .
```

... with [MP-Element][8] ...
```
docker build --build-arg MP_SRC=github_mpe -f docker/Dockerfile -t matpower/matpower:YYYY-MM-DD .
```

To build for a specific MATPOWER version, specify the corresponding GitHub tag
in `BRANCH` arg.
```
docker build --build-arg BRANCH=7.1 -f docker/Dockerfile -t matpower/matpower:7.1c .
```

To build for a specific Octave version, specify the tag in `BASE_TAG` arg.
```
docker build --build-arg BASE_TAG=5.2.0 --build-arg BRANCH=7.1 -f docker/Dockerfile -t matpower:7.1-oct-5.2.0 .
```

See `build_matpower_images.sh` for examples.

MATPOWER Legacy Images ([`matpower/matpower`][7]) (pre v7.0)
------------------------------------------------------------

See [`build_matpower_images.sh`] for examples.

---

[1]: https://hub.docker.com/r/matpower/octave
[2]: https://hub.docker.com/r/gnuoctave/octave
[3]: ./build_octave_images.sh
[4]: https://hub.docker.com
[5]: ./MATPOWER-Docker.md
[6]: ./Octave-Docker.md
[7]: https://hub.docker.com/r/matpower/matpower
[8]: https://github.com/MATPOWER/mp-element
