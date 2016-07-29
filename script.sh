#!/bin/bash

find . -type f | xargs sed -i 's/MyRawClusterBuilder/MyRawClusterBuilder/g'
