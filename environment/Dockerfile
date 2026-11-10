# Image with the pinned Python environment; the default command runs the quick checks.
# Build from the repository root:   docker build -f environment/Dockerfile -t second-order-ca .
# Run:                              docker run --rm second-order-ca
FROM python:3.12-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends gcc libc6-dev make bc \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /repo
COPY environment/requirements.txt environment/requirements.txt
RUN pip install --no-cache-dir -r environment/requirements.txt
COPY . .
CMD ["make", "check"]
