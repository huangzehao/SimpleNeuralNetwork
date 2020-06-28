# SimpleNeuralNetwork
This is a C++ implement of simple neural network. It's based on video [Neural Net in C++ Tutorial](https://vimeo.com/19569529) by David Miller.
# Test in Ubuntu
1 Gernerate training data to slove XOR problem
```
    g++ ./makeTrainingSamples.cpp -o makeTrainingSamples
    ./makeTrainingSamples > trainingData.
```
2 Test neural network
```
    g++ ./neural-net.cpp -o neural-net
    ./neural-net
```
And you will get the result!
# Tuning your model
3 Using different topology 
```
     In the trainingdata.txt, the default topology is 2 4 1.
     Changing to a different topology gives different results. 
     For example changing to 2 4 2 1 gives almost 30 percent more
     accurate model  for 4000 samples of data sets compared to previous one.
```
