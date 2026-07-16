FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN set -e \
      && apt-get -y update \
      && apt-get -y install --no-install-recommends --no-install-suggests \
      ca-certificates curl python3 samtools=1.16.1-1 python3-pip git

RUN set -e \
      && apt-get -y update \
      && apt-get install -y python3-setuptools \
      && pip3 install --no-cache-dir regex numpy==1.24.4 scipy==1.13.1 biopython==1.85 pandas==2.3.1 seaborn==0.13.2 matplotlib==3.10.5 \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*


RUN set -e \
      && cd / \
      && git clone https://github.com/nadeulab/IgCaller \
      && chmod +x IgCaller/IgCaller

ENV PATH=/IgCaller:${PATH}

ENTRYPOINT ["IgCaller"]
