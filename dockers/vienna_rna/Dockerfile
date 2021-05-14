FROM python:2.7

ENV VIENNA_DIR /usr/local

WORKDIR /root

# Update all packages
RUN apt-get update \
 && apt-get install -y build-essential

# ViennaRNA source
ADD https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz .

# Unpack and link ViennaRNA directory
RUN tar -xzvf ViennaRNA-2.4.17.tar.gz -C $VIENNA_DIR \
 && rm -f ViennaRNA-2.4.17.tar.gz

WORKDIR $VIENNA_DIR/ViennaRNA-2.4.17

# Configure and make ViennaRNA package
RUN ./configure \
 && make \
 && make install \
 && make clean

