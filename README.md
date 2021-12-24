# *Balanced-TLBO:*
#### *A Balanced Teaching Learning Based Optimization Algorithm*.
Teaching-learning-based optimization (TLBO) is a population-based metaheuristic algorithm which simulates the teaching and learning mechanisms in a classroom [1]. The TLBO algorithm has emerged as one of the most efficient and attractive optimization techniques. Even though the TLBO algorithm has an acceptable exploration capability and fast convergence speed, there may be a possibility to converge into a local optimum during solving complex optimization problems and there is a need to keep a balance between exploration and exploitation capabilities. 
Hence, a Balanced Teaching-Learning-Based Optimization (BTLBO) algorithm is proposed in this project. The proposed BTLBO algorithm is a modification of the TLBO algorithm and it consists of four phases: 
1) Teacher Phase in which a weighted mean is used instead of a mean value for keeping the diversity.
2) Learner Phase, which is same as the learner phase of basic TLBO algorithm.
3) Tutoring Phase, which is a powerful local search for exploiting the regions around the best ever found solution.
4) Restarting Phase, which improves exploration capability by replacing inactive learners with new randomly initialized learners.

## *CEC 2014 benchmark functions:*
Presented in a special session on real-parameter single objective optimization [2], the CEC 2014 benchmark test suite includes 30 optimization functions. These test functions can be categorized into four groups according to their characteristics.
- F01- F03 are uni-modal functions. This group has only one global optimum without any local optima and are suitable for examining the exploitation ability of algorithms.
- F04-F016 are simple multimodal functions and have many local optima. These are used to investigate the exploration and the local optima avoidance capabilities of algorithms.
- F17–F22 are hybrid functions and their subcomponents include both uni-modal and multi-modal functions. This group is suitable to examine exploitation and exploration of algorithms simultaneously.
- F23–F30 are composition functions which merge the features of the sub-functions better than the hybrid functions. These functions can maintain continuity around the local or global optima. The composition functions utilize the hybrid functions as the basic functions. In experiments, all these functions are rotated and shifted. Therefore, their complexity is increased dramatically.
