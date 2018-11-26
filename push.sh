#!/bin/bash
mvn clean
git add -A .
git commit -am "update"
git push iridium master
