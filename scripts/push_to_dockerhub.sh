#!/bin/bash

# Script to build and push cryvar Docker image to DockerHub
# Usage: ./push_to_dockerhub.sh [DOCKERHUB_USERNAME] [TAG]

set -e

DOCKERHUB_USERNAME=${1:-"dylanpilz"}
TAG=${2:-"latest"}
IMAGE_NAME="cryvar"

echo "Building Docker image..."
docker build -t ${IMAGE_NAME}:${TAG} .

echo "Tagging image for DockerHub..."
docker tag ${IMAGE_NAME}:${TAG} ${DOCKERHUB_USERNAME}/${IMAGE_NAME}:${TAG}

echo "Logging into DockerHub..."
docker login

echo "Pushing ${DOCKERHUB_USERNAME}/${IMAGE_NAME}:${TAG} to DockerHub..."
docker push ${DOCKERHUB_USERNAME}/${IMAGE_NAME}:${TAG}

echo "Successfully pushed ${DOCKERHUB_USERNAME}/${IMAGE_NAME}:${TAG} to DockerHub!"


