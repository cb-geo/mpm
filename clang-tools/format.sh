#!/bin/bash
find ../include/ -iname *.h -o -iname *.tcc | xargs clang-format -i
find ../src/ -iname *.h -o -iname *.cc | xargs clang-format -i
find ../tests/ -iname *.h -o -iname *.cc | xargs clang-format -i
