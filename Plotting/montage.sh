#!/bin/bash
DIR=graphics
montage -adjoin -tile 2x3 -geometry 864x432 \
  $DIR/$1_mean_enzyme_count_per_zone-Zone.1.svg $DIR/$1_mean_enzyme_count_per_zone-dxdt-Zone.1.svg \
  $DIR/$1_mean_enzyme_count_per_zone-Zone.2.svg $DIR/$1_mean_enzyme_count_per_zone-dxdt-Zone.2.svg \
  $DIR/$1_mean_enzyme_count_per_zone-Zone.3.svg $DIR/$1_mean_enzyme_count_per_zone-dxdt-Zone.3.svg \
  $DIR/$1-eiel.png

exit 0
