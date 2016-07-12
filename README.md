# Initial Performance Checks:
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


# Meeting notes:

## Ohio U:
* Have been trying to use past PHENIX clusterizer/cleaning it up. 
* Main goal is cluster splitting (?)
* Should see PHENIX clustering algorithm.
* Bremsstrahlung recovery step. 
    + Provide another clustering algorithm, works with CMS detector geometry.
* In CMS, five different clustering algorithms. 

## General Discussion.
* Learn about HIJING. 
    + embed single electron/photon events in HIJING events. 
    + simulate one shower in one particular place. 
    + ...efficient way of running without having to redo clustering all the time.
* Try sample of photons, sample of truth particles nearby. 
    + truth level
* Keep talking about obtaining samples, pi0 sample, etc.
* Something about ganging (?) the towers 2x2 and then you cluster. 

