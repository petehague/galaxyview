export PARENTDIR=/cosma5/data/dc-hagu1

rm -f datapath.txt
echo $PARENTDIR/survey_`date +%b%d%H%M` > datapath.txt
mkdir $PARENTDIR/survey_`date +%b%d%H%M`
bsub < image.job
