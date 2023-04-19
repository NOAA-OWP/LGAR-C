case=$1

if [ ! $# == 1 ]; then
    echo "Usage:  ./run_synthetic.sh OPTION={0,1,2}"
    echo "One of these options must be specified to run the corresponding synthetic example"
  exit
fi

if [ "${case}" == 0 ] || [ "${case}" == 1 ] || [ "${case}" == 2 ]; then
    echo 'Running synthetic example' ${case}
    ../build/lasam_standalone configs/config_lasam_synth_${case}.txt
else
    echo "Invalid option! ${case}"
    echo "Valid options: 0, 1, 2"
    exit
fi
    
