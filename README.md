# HASP
The implementation code of paper "A Spatio-temporal Optimized Framework for High-throughput
Shortest Distance Query Processing on Dynamic Road Networks" (Submitted to VLDB 2026). Please refer to the paper for the algorithm details. 


## Algorithms

The implementation code includes the index construction, query processing, and index update of our *DVPL* algorithm. 


## Data
The datasets of this paper are sourced from [NavInfo](https://en.navinfo.com/). We make the beijing dataset public in the `datasets` directory. After unzipping the zip file, you can run `./datasets beijing 6 1 32 NC 4 2 2 20 300 15 1 2 100 0.1 2` to have a taste.



## Dependency

1. `g++` and `boost`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
