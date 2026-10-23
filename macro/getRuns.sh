#!/bin/bash

#declare -a runNumbers=(79506 79508 79509 79510 79511 79512 79513 79514 79515 79516 79522 79523 79524 79525 79526 79527 79528 79529 79563 79564 79565 79566 79567 79568 79569 79570 79571 79572 80859 80860 80861 80864 80865 80866 80879 80880 80893 80895 80901 80902 80903 80904 80905 80906 80908 80909 80910 80912 80913 80915 81559 81560 81561 81562 81566 81576 81577 81580 81581 81583 81584 81585 81586 81600 81601 81602 81603 81605 81607 81608 81609 81610 81611 81612)
declare -a runNumbers=(79516)

build=ana578
tag=2026p004_v001
base=/sphenix/lustre01/sphnxpro/production/run3pp/physics/${build}_${tag}/DST_TRKR_TRACKS
for i in "${runNumbers[@]}"; do
  lo=$(( i / 100 * 100 )); hi=$(( lo + 100 ))
  dir=$base/$(printf 'run_%08d_%08d' $lo $hi)
  for j in {0..3}; do
    echo "Run number ${i} segment prefix ${j}"
    find "$dir" -name "DST_TRKR_TRACKS_run3pp_${build}_${tag}-$(printf '%08d' $i)-${j}*.root" | awk -F / '{print $NF}' | sort > jobLists_new/run${i}_${j}.list
  done
done
find . -maxdepth 1 -wholename 'jobLists_new/run*_*.list' -empty -delete   # drop empty lists

for i in "${runNumbers[@]}"; do
  cat jobLists_new/run${i}_*.list > jobLists_new/run${i}.list
  rm jobLists_new/run${i}_*.list
done
