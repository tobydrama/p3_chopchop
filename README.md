# CHOPCHOP
#### This repository is open sources as specified in the LICENSE file. It is AGPL-3.0.

#### Prerequisites:
- [Python](https://www.python.org/download/) 3.8
- [Docker](https://www.docker.com/products/docker-desktop/)
  - [Docker permissions](https://docs.docker.com/engine/install/linux-postinstall/). Several scoring methods are dockerized as they use python 2 specific libraries. To use them, the user / application must have permissions to run docker.
- Included libraries in [`lib/`](./lib/):
  - [Bowtie 1.0.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.1/)
  - [primer3](http://primer3.sourceforge.net/releases.php/)
  - [TwoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/)
  - ViennaRNA (docker)

#### Installation:
1. Clone the repository to your machine
2. Run `make` (takes a while).
3. Add paths to your genome directories to your `local_config.json` file.

#### Manual installation:
1. Install required python packages by running: `pip install -r requirements.txt`.
2. Make a copy of `config.json` & name it `config_local.json`.
3. Add paths to your genome directories to the config file.
