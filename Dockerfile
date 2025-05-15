FROM bentoml/model-server:0.11.0-py39
MAINTAINER ersilia

RUN pip install rdkit==2022.09.5
RUN pip install molecule-generation==0.4.1
RUN pip install tensorflow==2.9.1
RUN pip install numpy==1.23

WORKDIR /repo
COPY . /repo
