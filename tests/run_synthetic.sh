case=$1
echo 'Running synthetic example' ${case}
../build/lasam_standalone configs/config_lasam_synth_${case}.txt
