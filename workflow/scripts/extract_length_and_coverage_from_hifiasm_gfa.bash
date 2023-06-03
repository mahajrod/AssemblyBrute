#!/usr/bin/env bash

grep -P "^S" $1 | cut -f 2,4,5 | sed 's/rd:i://;s/LN:i://'
