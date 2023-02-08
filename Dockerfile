FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install molecule-generation
RUN conda install -c rdkit rdkit

WORKDIR /repo
COPY . /repo
