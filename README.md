# SimpleNeuralNetwork
This is a C++ implement of simple neural network. It's based on video [Neural Net in C++ Tutorial](https://vimeo.com/19569529) by David Miller.
# Test in Ubuntu
1 Gernerate training data to slove XOR problem
```
    g++ ./makeTrainingSamples.cpp -o makeTrainingSamples
    ./makeTrainingSamples > out.txt
```
2 Test neural netwrok
```
    g++ ./neural-net.cpp -o neural-net
    ./neural-net
```
And you will get the result!

