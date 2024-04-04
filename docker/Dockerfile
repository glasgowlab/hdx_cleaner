FROM jupyter/scipy-notebook

COPY --chown=jovyan . /app/PIGEON_FEATHER
RUN cd /app/PIGEON_FEATHER && pip install .
RUN cd /app/ && git clone https://github.com/lucl13/bayesian_hdx.git 

RUN pip install --index-url https://pypi.cs.uni-tuebingen.de/simple/ pyopenms
RUN pip install HDXrate
RUN conda install conda-forge::pymol-open-source

# install fonts
USER root
COPY fonts /usr/share/fonts/truetype/myfonts
RUN apt-get update && apt-get install -y fontconfig && fc-cache -f -v
# refresh matplotlib font cache
RUN rm -fr ~/.cache/matplotlib

USER jovyan
WORKDIR /work