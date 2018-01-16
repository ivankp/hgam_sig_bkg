#!/bin/bash

lumi=36067
window=0.170212766

for b in \
  'pT_yy:100:0,500' \
  'pT_yy:0,10,20,30,45,60,80,120,170,200,350,500'
do

coarse=$(echo $b | sed -r 's/[^:,]+:[^:,]+:[^:]+//' | sed 's/.\+/_coarse/')
pref=hist$coarse

./bin/hist \
  -i ../flat/hgam_data_h015d.root ../flat/hgam_mc_h015d_* \
  -o ${pref}.root -b $b

rxplot ${pref}.root -o ${pref}_split.pdf --logy \
  -r 'sid/all_MC' \
     'md/.*ybyt/& (#minus )/scale=-1,' \
     "md/(.*)data(.*)/\\1background (data)\\2/scale=${window}," \
     "md/h\\d+\\l?_.*//scale=${lumi}," \
     'dt/.*(h\d+\l?).*/\1' 'tg' 'n+t/^.*/& ' \
     'dl/h\d+\l?/' 'll/^_/' 'll/_$/' 'l/_/ ' \
     'nx/^.*/& [GeV]' \
     'y/^.*/Num events' \
   -ltr \
   --colors=602 46 8 90 93 95 6 13 12 28

rxplot ${pref}.root -o ${pref}_sum.pdf --logy \
  -r 'sd/.*(all_MC|data).*' \
     "md/.*data.*/background (data)/scale=${window}," \
     "md/.*_MC.*/signal (MC)/scale=${lumi}," \
     'dt/.*(h\d+\l?).*/\1' 'n+t/^.*/& ' \
     'gg/^.*/' \
     'dl/h\d+\l?/' 'll/^_/' 'll/_$/' 'l/_/ ' \
     'nx/^.*/& [GeV]' \
     'y/^.*/Num events' \
   -ltr \
   --colors=602 46 \
   ${coarse:+--val-fmt=.0f}

done
