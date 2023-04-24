python3 gen-input.py -n 30000 -m 4 -z 230 -o input1
python3 gen-input.py -n 30000 -m 4 -z 330 -o input2
python3 gen-output.py -f input1 input2 -e 2 -o output1
nvcc main2.cu
#g++ main.cpp
./a.out input1 input2 output2
#python3 gen-output.py -f input2 input1 -e 2 -o output2
python3 checker.py -f output1 output2 -e 4
#python3 checker.py -f output1 output2 -e 4
