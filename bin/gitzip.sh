#!/bin/bash
GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
GIT_VERSION=`git describe --abbrev=16 --dirty --always`
FILENAME=$GIT_BRANCH-$GIT_VERSION.zip
echo "Creating archive $FILENAME"
git archive -o $FILENAME HEAD
