FROM node:22-bookworm AS frontend-builder
WORKDIR /app/frontend
COPY frontend/package*.json ./
RUN npm ci
COPY frontend/ ./
RUN npm run build

FROM golang:1.26.10-bookworm AS backend-builder
WORKDIR /app/src
COPY src/go.mod src/go.sum ./
RUN go mod download
COPY src/main.go ./
RUN CGO_ENABLED=0 go build -trimpath -o /out/rmsdapp main.go

FROM debian:bookworm-slim AS estimator-builder
RUN apt-get update && apt-get install -y --no-install-recommends g++ libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /app/src
COPY src/estimate_hinge_numbers.cpp src/*.h ./
RUN g++ estimate_hinge_numbers.cpp -o /out-estimator -std=c++14 -Wall -Wextra \
    -O3 -Werror -fopenmp -fstack-protector-strong -D_FORTIFY_SOURCE=2 -Wl,-z,relro,-z,now

# Optional research environment; excluded from the production image.
FROM python:3.12-slim-bookworm AS experiments
RUN apt-get update && apt-get install -y --no-install-recommends g++ libeigen3-dev git \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir poetry==2.2.1
WORKDIR /app/src
COPY src/ ./
CMD ["bash"]

FROM debian:bookworm-slim AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates libgomp1 libstdc++6 \
    && rm -rf /var/lib/apt/lists/* \
    && useradd --create-home --uid 10001 --shell /usr/sbin/nologin app
WORKDIR /app/src
COPY --from=backend-builder /out/rmsdapp ./rmsdapp
COPY --from=estimator-builder /out-estimator ./estimate_hinge_numbers
COPY src/template.html ./
COPY --from=frontend-builder /app/frontend/dist /app/frontend/dist
ENV GIN_MODE=release
USER app
EXPOSE 8080
CMD ["./rmsdapp"]
