FROM python:3.12-slim

# Add bash for WDL/Cromwell entry scripts and OpenMP lib for scikit-learn
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash \
      libgomp1 \
 && rm -rf /var/lib/apt/lists/*

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Put everything under /build/PhylogicNDT
ARG APP_DIR=/build/PhylogicNDT
WORKDIR ${APP_DIR}

# 1) Install deps first for layer caching
COPY requirements.txt /tmp/requirements.txt
RUN python -m pip install --upgrade pip wheel setuptools \
 && python -m pip install --no-cache-dir --prefer-binary -r /tmp/requirements.txt

# 2) Copy your repo and install your package
COPY . ${APP_DIR}
RUN python -m pip install --no-cache-dir .
CMD ["/bin/bash"]
