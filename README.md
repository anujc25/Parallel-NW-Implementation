# Parallel-NW-Implementation
Parallel implementation of NW algorithms with NVIDIA GPU and CUDA C++


This work is focused on optimization of well known DNA Sequence Assembly tool Velvet and DNA Alignment algorithm Needleman-Wunsch using CUDA and HPC technologies. This project focus on performance comparison of different implementation of NW algorithm.

This is new parallel approach of Needleman-Wunsch algorithm for global sequence alignment. This approach uses skewing transformation for traversal and calculation of the dynamic programming matrix. We compare the execution time of sequential CPU based implementation with two parallel GPU based implementations: Single-kernel invocation with lock-free block synchronization and multi-kernel invocation at block-synchronization points. Both the GPU based implementations gave upto 6 times performance improvement over the sequential CPU based implementation.

For more information refer to the following paper on "A GPU based implementation of Needleman-Wunsch algorithm using skewing transformation"
Link: http://ieeexplore.ieee.org/document/7346733/
