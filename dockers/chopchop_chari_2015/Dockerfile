FROM python:2.7

WORKDIR /chopchop_chari_2015

RUN pip install --upgrade pip \
 && pip install scipy


ADD . .

ENTRYPOINT ["python", "./score_chari_2015.py"]