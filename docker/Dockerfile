FROM jupyter/scipy-notebook

RUN conda install jupyterlab mdanalysis numba
RUN conda install pymol-open-source
RUN pip install pyopenms hdxrate 

COPY --chown=jovyan . /app/PIGEON_FEATHER
RUN cd /app/PIGEON_FEATHER && pip install .
RUN cd /app/ && git clone https://github.com/lucl13/bayesian_hdx.git && cd bayesian_hdx && pip install .

# install fonts
USER root
COPY fonts /usr/share/fonts/truetype/myfonts
RUN apt-get update && apt-get install -y fontconfig && fc-cache -f -v
# refresh matplotlib font cache
RUN rm -fr ~/.cache/matplotlib

USER jovyan
WORKDIR /work