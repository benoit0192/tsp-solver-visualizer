DOCKER_DIR=$(dirname "$0")
IMAGE_TAG=tsp:latest
CONTEXT_DIR=$DOCKER_DIR/.

docker build                        \
    -f $DOCKER_DIR/tsp.Dockerfile   \
    -t $IMAGE_TAG                   \
    $CONTEXT_DIR
