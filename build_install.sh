#!/bin/bash
R CMD build .
sudo R COMD REMOVE --library=/usr/libg/R/library/ rpart
sudo R CMD INSTALL --library=/usr/lib/R/library/ $1 
