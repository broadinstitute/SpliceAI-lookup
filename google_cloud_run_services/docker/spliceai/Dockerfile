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

RUN python3 -m pip install tensorflow==2.16.1

RUN apt update && apt-get install --no-install-recommends -y build-essential libpq-dev

COPY docker/spliceai/requirements.txt /
RUN python3 -m pip install --upgrade -r /requirements.txt

ARG RANDOM=2
RUN python3 -m pip install https://github.com/bw2/SpliceAI/archive/refs/heads/master.zip

ARG CONCURRENCY="2"
ARG GENOME_VERSION="unknown"

COPY docker/ref/GRCh${GENOME_VERSION} /
COPY docker/spliceai/annotations/GRCh${GENOME_VERSION} /
COPY server.py /

ENV TF_CPP_MIN_LOG_LEVEL=3

ENV PORT=8080
ENV TOOL=spliceai
ENV GENOME_VERSION=${GENOME_VERSION}
ENV CONCURRENCY=${CONCURRENCY}
ENV RUNNING_ON_GOOGLE_CLOUD_RUN=1

CMD exec gunicorn --preload --bind :$PORT --workers ${CONCURRENCY} --threads 1 --timeout 0 server:app
