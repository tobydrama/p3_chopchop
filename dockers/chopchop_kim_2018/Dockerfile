FROM python:2.7

WORKDIR /chopchop_kim_2018

RUN pip install --upgrade pip \
 && pip install keras==2.3 np_utils theano enum34

ADD . .

ENTRYPOINT ["python", "./score_kim_2018.py"]


