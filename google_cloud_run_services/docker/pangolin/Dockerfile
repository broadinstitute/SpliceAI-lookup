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
#
# --extra-index-url, not --index-url: the latter REPLACES PyPI, so any build
# dependency pip needs that the PyTorch index doesn't carry is simply
# unfindable ("No matching distribution found for flit_core"). That failure
# stayed hidden for as long as this layer was restored from the build cache
# and only appeared once a cache miss forced a cold rebuild. Keeping PyPI in
# the search path fixes it, and the explicit "+cpu" local version still pins
# the CPU wheel unambiguously, since only the PyTorch index publishes it.
RUN python3 -m pip install torch==2.2.1+cpu --extra-index-url https://download.pytorch.org/whl/cpu

RUN apt update && apt-get install --no-install-recommends -y build-essential libpq-dev

COPY docker/pangolin/requirements.txt /
RUN python3 -m pip install --upgrade -r /requirements.txt

# Pinned to a commit rather than whatever the default branch happens to be at build time, so a
# rebuild cannot silently change the scores this service returns. See the equivalent comment in
# docker/spliceai/Dockerfile. Bump deliberately, then re-run test_api_consistency.py against
# the dev deployment.
ARG PANGOLIN_COMMIT=31b5aef6ad74a09b630443d9672b830ec4947bb6
RUN git clone https://github.com/bw2/Pangolin.git \
    && cd Pangolin \
	&& git checkout ${PANGOLIN_COMMIT} \
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

# --timeout 870 (was 290, 120, 0): recycle a worker stuck >870s. Raised 290 -> 870 to match
# spliceai (with Cloud Run's request timeout raised 300 -> 900 in build_and_deploy.py, so
# gunicorn recycles first for a clean 500 rather than a platform 504) so heavy inferences
# complete instead of being killed. Trade-off: a sync worker (1 of WORKERS) is tied up for
# up to 870s on such a request, lowering the instance's effective concurrency while it runs.
# --graceful-timeout gives in-flight requests 30s to finish on restart.
#
# No --preload: kept identical to the spliceai image so both tools load their models in the
# worker rather than in the arbiter before the fork. See that Dockerfile and preload_models().
CMD exec gunicorn --bind :$PORT --workers ${WORKERS} --threads 1 --timeout 870 --graceful-timeout 30 server:app
