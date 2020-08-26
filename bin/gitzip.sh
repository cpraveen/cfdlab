#!/bin/bash
GIT_BRANCH=`git rev-parse --abbrev-ref HEAD`
GIT_VERSION=`git describe --abbrev=16 --dirty --always`
BASEDIR=`basename "$PWD"`
FILENAME=$BASEDIR-$GIT_BRANCH-$GIT_VERSION.zip
echo "Creating archive $FILENAME"
git archive --prefix=$BASEDIR/ -o $FILENAME HEAD
