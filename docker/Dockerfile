FROM ubuntu
LABEL Chenlin Lu <lucltsinghua@gmail.com>

RUN apt-get update && apt-get install -y git 

FROM jupyter/scipy-notebook
RUN pip install --index-url https://pypi.cs.uni-tuebingen.de/simple/ pyopenms
RUN pip install HDXrate

RUN git clone https://github.com/glasgowlab/PIGEON.git
RUN git clone https://github.com/lucl13/bayesian_hdx.git 
WORKDIR /home/jovyan/work