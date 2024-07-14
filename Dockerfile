FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y sudo && \
    apt-get install -y build-essential && \
    apt-get install -y emacs && \
    apt-get install -y wget && \
    apt-get install -y zlib1g && \
    apt-get install -y unzip && \
    apt-get install -y git && \
    apt-get install -y cmake && \
    apt-get install -y libssl-dev && \
    apt-get install -y libz-dev
RUN apt-get install -y gcc && \
    apt-get install -y g++
ENV CC=/usr/bin/gcc \
    CXX=/usr/bin/g++
RUN apt-get install -y libeigen3-dev
RUN wget https://dl.google.com/go/go1.22.5.linux-amd64.tar.gz && \
    tar -C /usr/local -xzf go1.22.5.linux-amd64.tar.gz && \
    rm go1.22.5.linux-amd64.tar.gz
ENV PATH=$PATH:/usr/local/go/bin \
    GOPATH=$HOME/go
WORKDIR $HOME/usr/src
COPY src/go.mod .
COPY src/go.sum .
COPY src/estimate_hinge_numbers.cpp .
COPY src/*.h .
RUN go mod download
RUN g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -mtune=native -march=native -mfpmath=both -Werror -fopenmp
RUN apt-get update && \
    apt-get install -y python3 python3-pip
RUN wget -O- https://install.python-poetry.org | python3 -
ENV PATH=$PATH:/root/.local/bin

CMD bash
