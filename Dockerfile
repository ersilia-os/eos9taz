FROM bentoml/model-server:0.11.0-py310
MAINTAINER ersilia

RUN pip install molecule-generation==0.4.1
RUN pip install numpy==1.26.4
RUN pip install rdkit==2025.9.1

WORKDIR /repo
COPY . /repo
