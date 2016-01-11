#!/bin/bash

# Script to download the latest commit from GitHub
# Assumes that you have SSH keys set up. 
# Run in directory with subdirectories as repos named after SHA. 

GIT_URL="$1"

git clone "$GIT_URL" clone
cd clone	

LATEST_TAG=$(git describe --abbrev=0 --tags)
echo "Latest tag is $LATEST_TAG"

if [ ! -d "../$LATEST_TAG" ]; then
	echo "Keeping latest tagged version..."
	git checkout -b "tags/$LATEST_TAG" "tags/$LATEST_TAG"
	cd ..
	echo "Renaming repo..."
	mv clone "$LATEST_TAG"
	echo "Updating `latest` symlink..."
	rm "latest"
	ln -s "$LATEST_TAG" "latest"
else
	echo "Latest tagged version already present..."
	echo "Deleting duplicate repo..."
	cd ..
	rm -rf clone
fi
