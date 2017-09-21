#!/bin/bash
set -v

pip install git+https://github.com/${TRAVIS_REPO_SLUG%/*}/cobrame.git@${TRAVIS_BRANCH}