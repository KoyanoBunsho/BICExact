FROM node:22-bookworm AS frontend-builder

WORKDIR /app/frontend
COPY frontend/package*.json ./
RUN npm install
COPY frontend/ ./
RUN npm run build

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
    apt-get install -y libz-dev && \
    apt-get install -y golang-go
RUN apt-get install -y gcc && \
    apt-get install -y g++
ENV CC=/usr/bin/gcc \
    CXX=/usr/bin/g++
RUN apt-get install -y libeigen3-dev
ENV GOPATH=/root/go
WORKDIR /root/usr/src
COPY src/go.mod .
COPY src/go.sum .
COPY src/estimate_hinge_numbers.cpp .
COPY src/*.h .
COPY src/main.go .
COPY src/template.html .
COPY --from=frontend-builder /app/frontend/dist /root/usr/frontend/dist
RUN go mod download
RUN g++ estimate_hinge_numbers.cpp -o estimate_hinge_numbers -std=c++14 -lstdc++fs -Wall -Wextra -O3 -Werror -fopenmp
RUN go build -o rmsdapp main.go
RUN apt-get update && \
    apt-get install -y python3 python3-pip
RUN wget -O- https://install.python-poetry.org | python3 -
ENV PATH=$PATH:/root/.local/bin

EXPOSE 8080
CMD ["./rmsdapp"]
