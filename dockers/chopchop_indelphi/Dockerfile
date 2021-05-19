FROM python:2.7

WORKDIR /chopchop_inDelphi

RUN pip install --upgrade pip \
 && pip install scikit-learn==0.18.1 numpy pandas scipy

# Add script
ADD . .

ENTRYPOINT ["python", "./run_repair_prediction.py"]