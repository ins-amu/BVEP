FROM centos:7

ENV stan_ver=2.28.2

RUN yum install centos-release-scl -y
RUN yum install devtoolset-10 -y

WORKDIR /work

RUN curl -LO https://github.com/stan-dev/cmdstan/releases/download/v${stan_ver}/cmdstan-${stan_ver}.tar.gz \
 && tar xzf cmdstan-${stan_ver}.tar.gz \
 && rm $fname \
 && cd cmdstan-$stan_ver \
 && scl enable devtoolset-10 "make build -j1"

