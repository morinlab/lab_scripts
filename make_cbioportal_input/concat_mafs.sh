#!/bin/bash

(>&2 echo "Concatenating all MAFs...")

egrep -h "^(#|Hugo_Symbol)" $1;
egrep -h -v "^(#|Hugo_Symbol)" $@
