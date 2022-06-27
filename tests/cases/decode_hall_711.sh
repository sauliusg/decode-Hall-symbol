#!/bin/sh

./bin/decode_hall --version | awk '{print $1, $2}'
