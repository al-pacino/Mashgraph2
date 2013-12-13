mkdir -p build
cd build
echo --------------------------------------------------------------------------
cmake -G "Unix MakeFiles" ../
echo --------------------------------------------------------------------------
make
echo --------------------------------------------------------------------------
task2 -d ../data/train_labels.txt -m model.txt --train
echo --------------------------------------------------------------------------
task2 -d ../data/test_labels.txt -m model.txt -l predictions.txt --predict
echo --------------------------------------------------------------------------
echo
task2_test -g ../data/test_labels.txt -p predictions.txt
