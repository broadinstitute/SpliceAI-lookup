FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    wget \
    bzip2 \
    unzip \
    git \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev

RUN python3 -m pip install torch==2.2.1+cpu -f https://download.pytorch.org/whl/torch_stable.html

RUN apt update && apt-get install --no-install-recommends -y build-essential libpq-dev

COPY docker/pangolin/requirements.txt /
RUN python3 -m pip install setuptools==57.5.0 \
    && python3 -m pip install pyvcf>=0.6.8 \
	&& python3 -m pip install --upgrade setuptools \
	&& python3 -m pip install --upgrade -r /requirements.txt

RUN git clone https://github.com/bw2/Pangolin.git \
    && cd Pangolin \
	&& python3 -m pip install .

ARG GENOME_VERSION="unknown"

COPY docker/ref/GRCh${GENOME_VERSION} /
COPY docker/pangolin/annotations/GRCh${GENOME_VERSION} /
COPY server.py /

ENV PORT=8080
ENV TOOL=pangolin
ENV GENOME_VERSION=${GENOME_VERSION}

CMD exec gunicorn --preload --bind :$PORT --workers 1 --threads 1 --timeout 0 server:app
