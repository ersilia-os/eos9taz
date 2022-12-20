FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit=2021.03.4
RUN pip install joblib==1.1.0
RUN pip install tensorflow==2.9.1
RUN pip install dpu-utils == 0.2.13
RUN pip install more-itertools
RUN pip install numpy
RUN pip install pandas
RUN pip install protobuf==3.19.5
RUN pip install scikit-learn==0.24.1
RUN pip install tf2_gnn==2.13.0
RUN pip install SQLAlchemy==1.3.24


WORKDIR /repo
COPY . /repo
