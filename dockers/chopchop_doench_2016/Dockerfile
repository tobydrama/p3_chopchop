FROM python:2.7

WORKDIR /chopchop_doench_2016

# biopython==1.76 is the last release to work with python2.7,
# sckit-learn==0.18.1 is required for the doench_2016 model.
RUN pip install --upgrade pip \
 && pip install pandas numpy scipy biopython==1.76 scikit-learn==0.18.1

# Add doench_2016 models to the WORKDIR
ADD models/ models/

# Add script
ADD . .

ENTRYPOINT ["python", "./score_doench_2016.py"]