#!/bin/bash
for file in * ; do lower=$(echo $file | tr A-Z a-z) && [[ $lower != $file ]] && mv $file $lower ; done
