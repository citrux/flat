#!/bin/bash

find . -regex ".*\(cc\|hh\)" | xargs clang-format -style=file -i
