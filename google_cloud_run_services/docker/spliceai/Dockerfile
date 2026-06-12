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

# tensorflow-cpu (vs the default tensorflow wheel) skips the bundled CUDA
# libraries (~1.5 GB) we never use on Cloud Run — `import tensorflow` and the
# rest of the API are unchanged.
RUN python3 -m pip install tensorflow-cpu==2.16.1

RUN apt update && apt-get install --no-install-recommends -y build-essential libpq-dev

COPY docker/spliceai/requirements.txt /
RUN python3 -m pip install --upgrade -r /requirements.txt

ARG RANDOM=2
RUN python3 -m pip install https://github.com/bw2/SpliceAI/archive/refs/heads/master.zip

ARG WORKERS="2"
ARG GENOME_VERSION="unknown"

COPY docker/ref/GRCh${GENOME_VERSION} /
COPY docker/spliceai/annotations/GRCh${GENOME_VERSION} /
COPY server.py /
COPY sai10k_predictions.py /

ENV TF_CPP_MIN_LOG_LEVEL=3

# Cap every CPU thread pool to 1 thread per worker. Without these, each gunicorn
# worker's TensorFlow/Eigen/OpenMP/BLAS pools default to the instance core count, so
# WORKERS workers on a 2-CPU instance oversubscribe the CPUs by a large factor. Under a
# burst of cache-miss requests a single ~6s inference then balloons into the minutes and
# trips the 300s Cloud Run timeout (504). One thread per worker keeps per-request CPU
# bounded; parallelism comes from running WORKERS workers, not threads within a worker.
ENV OMP_NUM_THREADS=1
ENV TF_NUM_INTRAOP_THREADS=1
ENV TF_NUM_INTEROP_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

ENV PORT=8080
ENV TOOL=spliceai
ENV GENOME_VERSION=${GENOME_VERSION}
ENV WORKERS=${WORKERS}
ENV RUNNING_ON_GOOGLE_CLOUD_RUN=1

# --timeout 870 (was 290, 120, 0): recycle a worker stuck >870s. Raised 290 -> 870 — with
# Cloud Run's request timeout raised 300 -> 900 in build_and_deploy.py, so gunicorn still
# recycles first (clean 500) rather than the platform (504) — so the heaviest distance=10000
# requests on many-isoform genes, which can exceed 290s under CPU contention, complete
# instead of being killed. Trade-off: a sync worker (1 of WORKERS) is tied up for up to 870s
# on such a request, lowering the instance's effective concurrency while it runs.
# --graceful-timeout gives in-flight requests 30s to finish on restart.
CMD exec gunicorn --preload --bind :$PORT --workers ${WORKERS} --threads 1 --timeout 870 --graceful-timeout 30 server:app
