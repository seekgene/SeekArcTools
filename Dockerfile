FROM ubuntu:22.04

# Prevent interactive prompts during apt installations
ENV DEBIAN_FRONTEND=noninteractive

# Use Aliyun's official mirror for Ubuntu apt-get
RUN sed -i 's@//.*archive.ubuntu.com@//mirrors.aliyun.com@g' /etc/apt/sources.list && \
    sed -i 's@//.*security.ubuntu.com@//mirrors.aliyun.com@g' /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends curl bzip2 ca-certificates git build-essential && \
    rm -rf /var/lib/apt/lists/*

# Install micromamba manually using robust method
RUN curl -fSL --retry 3 --connect-timeout 10 --max-time 120 \
        https://micro.mamba.pm/api/micromamba/linux-64/latest \
        -o /tmp/micromamba.tar.bz2 && \
    tar -xvjf /tmp/micromamba.tar.bz2 -C /usr/local/ bin/micromamba && \
    rm -f /tmp/micromamba.tar.bz2 && \
    echo "micromamba installed successfully at /usr/local/bin/micromamba"

# Initialize micromamba environment variables
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:/usr/local/bin/micromamba:$PATH

# Set the working directory
WORKDIR /app

# Copy the environment file first to leverage Docker cache
COPY env.yaml /app/env.yaml

# Configure pip to use Aliyun mirror
ENV PIP_INDEX_URL=https://mirrors.aliyun.com/pypi/simple/

# Install the dependencies using micromamba
RUN micromamba install -y -n base -f /app/env.yaml && \
    micromamba clean --all -y

# Copy the rest of the application code
COPY . /app

# Install the Python package itself
RUN pip install --no-cache-dir .

# Provide a fallback wrapper just in case scripts call /seekarctools/seekarctools
RUN mkdir -p /seekarctools && ln -s /opt/conda/bin/seekarctools /seekarctools/seekarctools

# Define the entrypoint
ENTRYPOINT ["seekarctools"]
