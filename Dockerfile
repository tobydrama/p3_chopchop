#FROM ubuntu:latest
#WORKDIR /root
#RUN apt-get update && apt-get install -y coreutils gcc wget apt-utils libgsl23 && \
#    wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_20_04/viennarna_2.4.17-1_amd64.deb && \
#    dpkg -i viennarna_2.4.17-1_amd64.deb

FROM ubuntu:latest

RUN apt-get update && \
	apt-get install -y apt-utils wget g++ libmpfr-dev

WORKDIR /viennaRNA

RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.17.tar.gz && \
	tar -zxvf ViennaRNA-2.4.17.tar.gz && \
	ViennaRNA-2.4.17/configure
	# && \
#	make && make install && \
#	rm ViennaRNA-2.4.17.tar.gz

FROM python:2.7

COPY CRISPR_Calculations.py .
COPY CRISPRoff/ ./CRISPRoff/


#ENTRYPOINT ["python","./CRISPR_Calculations.py"]
