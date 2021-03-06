# General Overview:
This an sPHENIX module that implements the island algorithm (the 'IA' in 'RawClusterBuilderIA') as used by the CMS
collaboration for clustering in the EMCal. The algorithm itself can be found in include/IslandAlgorithm.h.

# Bugs discovered:
* [__FIXED__] Forgot to ensure that no double-counting of towers occurs. In other words, I found events containing two
clusters that have one or more towers in common. 

# Initial Performance Checks (old):
(threshold seed energy set at 0.1 GeV). 

## __Single Electron__:
|                       | Input | Output 1  | Output 2 | Output 3  |
| --------------------- | ----- | --------- | -------- | --------  |
| __Energy__            | 20.   | 18.05     | 19.26    | 19.42     |
| __Eta__               | 0.0   | 0.04      | 0.01     | -0.01     |
| __Phi__               | 0.0   | -0.01     | -0.01    | -0.01     |

## __Single Photon__:
|                       | Input | Output 1  | Output 2  | Output 3  |
| --------------------- | ----- | --------- | --------  | --------  |
| __Energy__            | 20.   | 17.39     | 15.54     | 17.39     |
| __Eta__               | 0.0   | -0.04     | -0.01     | -0.03     |
| __Phi__               | 0.0   | 0.00      | 0.00      | 0.00     |


## __Single Pi0__:
|                       | Input | Output 1  | Output 2  | Output 3  |
| --------------------- | ----- | --------- | --------  | --------  |
| __Energy__            | 20.   | 11.94     | 16.77     | 15.98     |
| __Eta__               | 0.0   | 0.01     | 0.03     | 0.03     |
| __Phi__               | 0.0   | 0.00     | 0.00      | 0.00     |

