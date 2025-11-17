#! /usr/bin/env bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t bguo068/ishare:v0.1.11 \
  -t bguo068/ishare:latest \
  --push .
