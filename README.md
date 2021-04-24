# CHOPCHOP
#### This repository is open sources as specified in the LICENSE file. It is AGPL-3.0.

#### Prerequisites:
- [Python](https://www.python.org/download/) 3.8
- [Docker](https://www.docker.com/products/docker-desktop/)
-- [Docker permissions](https://docs.docker.com/engine/install/linux-postinstall/). Scoring methods Alkan 2018 & Doench 2016 are dockerized as they use python 2 specific libraries. To use them, the user must have permissions to run docker.
- Included libraries in [`lib/`](./lib/):
-- [Bowtie 1.0.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.1/)
-- [primer3](http://primer3.sourceforge.net/releases.php/)
-- [TwoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/)

#### Installation:
1. Clone the repository to your machine
2. Install required python packages by running: `pip install -r requirements.txt`.
2. Make a copy of `config.json` & name it `config_local.json`.
3. Add paths to your genome directories to the config file.
4. If you plan on using Alkan 2018 & Doench 2016 scoring methods, [build the docker containers](#docker-containers)

##### Docker containers:
To use scoring methods Alkan 2018 & Doench 2016, their respective docker containers must be built by running:
```sh
cd /path/to/chopchop/dockers/chopchop_crisproff  # Alkan 2018
docker build -t chopchop_cripsroff .  # This takes a very long time
```

```sh
cd /path/to/chopchop/dockers/chopchop_doench_2016  # Doench 2016
docker build -t chopchop_doench_2016 .
```
