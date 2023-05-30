FROM ruby:2.7.0

LABEL maintainer="BarKeeper (Kai Müller, Sarah Wiechers)"

RUN addgroup --gid 1000 barkeeper
RUN adduser --disabled-password --gecos '' --uid 1000 --gid 1000 barkeeper

ARG PUMA_PORT
ARG RAILS_ENV

RUN apt-get update -qq && apt-get install -y build-essential libpq-dev nodejs cmake

### Setup for BarPipe

# Install TRE library
RUN apt-get install -y tre-agrep libtre5 libtre-dev

# Install BLAST
RUN apt-get install -y ncbi-blast+

# Copy Usearch to container
COPY lib/usearch /usr/local/bin/usearch
ENV PATH "$PATH:/usr/local/bin/usearch"

# Conda setup

# Install base utilities
RUN apt-get update \
    && apt-get install -y wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install Conda packages used in BarPipe
RUN conda install -c bioconda lima
RUN conda install -c bioconda pbccs
RUN conda install -c bioconda samtools
RUN conda install -c bioconda seqtk

###

ENV RAILS_ROOT /var/www/barkeeper
RUN mkdir -p $RAILS_ROOT
WORKDIR $RAILS_ROOT

COPY Gemfile Gemfile
COPY Gemfile.lock Gemfile.lock

ENV BUNDLER_VERSION=2.3.5
RUN gem install rails bundler:2.3.5

RUN if [ "$RAILS_ENV" = "development" ]; then bundle config set --local without test; else bundle config set --local without test:development; fi
RUN bundle install
RUN chown -R barkeeper:barkeeper $RAILS_ROOT

RUN rails assets:precompile
EXPOSE $PUMA_PORT

CMD ["bundle", "exec", "puma", "-C", "config/puma.rb"]
