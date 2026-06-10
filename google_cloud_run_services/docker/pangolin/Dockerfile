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

# Pull torch from PyTorch's CPU-only wheel index — the default PyPI wheel
# bundles CUDA/cuDNN/cuBLAS (~2 GB on disk) that we never use on Cloud Run.
# Important: this must come before the requirements.txt install below; if
# pip first resolves torch from PyPI it will pick the GPU build.
RUN python3 -m pip install torch==2.2.1 --index-url https://download.pytorch.org/whl/cpu

RUN apt update && apt-get install --no-install-recommends -y build-essential libpq-dev

COPY docker/pangolin/requirements.txt /
RUN python3 -m pip install --upgrade -r /requirements.txt

RUN git clone https://github.com/bw2/Pangolin.git \
    && cd Pangolin \
	&& python3 -m pip install .

ARG WORKERS="2"
ARG GENOME_VERSION="unknown"

COPY docker/ref/GRCh${GENOME_VERSION} /
COPY docker/pangolin/annotations/GRCh${GENOME_VERSION} /
COPY server.py /
COPY sai10k_predictions.py /

# Cap every CPU thread pool to 1 thread per worker. Without these, each gunicorn
# worker's torch/OpenMP/BLAS pools default to the instance core count, so WORKERS workers
# on a 2-CPU instance oversubscribe the CPUs by a large factor and a burst of cache-miss
# inferences balloons past the 300s Cloud Run timeout (504). One thread per worker keeps
# per-request CPU bounded; parallelism comes from WORKERS workers, not threads per worker.
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

ENV PORT=8080
ENV TOOL=pangolin
ENV GENOME_VERSION=${GENOME_VERSION}
ENV WORKERS=${WORKERS}
ENV RUNNING_ON_GOOGLE_CLOUD_RUN=1

# --timeout 290 (was 120, was 0): recycle a worker stuck >290s so a saturated instance
# self-heals — `--timeout 0` never recycled a wedged worker, so 504s piled onto one
# instance for hours. Raised 120 -> 290 to match spliceai, so heavy inferences finish
# instead of being killed at 120s. Kept just under Cloud Run's 300s request timeout so
# gunicorn recycles the worker first (clean 500) rather than the platform (504).
# --graceful-timeout gives in-flight requests 30s to finish on restart.
CMD exec gunicorn --preload --bind :$PORT --workers ${WORKERS} --threads 1 --timeout 290 --graceful-timeout 30 server:app
