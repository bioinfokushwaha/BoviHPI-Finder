FROM ubuntu:22.04

# Install micromamba
RUN apt-get update && apt-get install -y \
    wget bzip2 ca-certificates curl git libglib2.0-0 libxext6 libsm6 libxrender1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    MAMBA_EXE=/bin/micromamba \
    PATH=/opt/conda/bin:$PATH

RUN curl -L https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C /bin/ --strip-components=1 bin/micromamba

# Create env and install dependencies via micromamba
COPY environment.yml .
RUN micromamba create -y -n bovihpi -f environment.yml && \
    micromamba clean --all --yes

# Set environment
SHELL ["micromamba", "run", "-n", "bovihpi", "/bin/bash", "-c"]

# Copy pipeline scripts and tools
WORKDIR /opt/project
COPY Exe/ Exe/
COPY Scripts/ Scripts/
COPY master.sh .

# Make main script executable
RUN chmod +x master.sh

CMD ["micromamba", "run", "-n", "bovihpi", "--"]
CMD ["./master.sh"]

