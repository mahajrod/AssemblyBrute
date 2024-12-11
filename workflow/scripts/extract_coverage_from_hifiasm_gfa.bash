#!/usr/bin/env bash

grep -P "^S" $1 | cut -f 2,5 | sed 's/rd:i://' | awk '{print $1"\t"$2+1}'
