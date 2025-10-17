# vesicle-ml-shapes

Directory Structure
```
project_root/
├─ bvp6c-solver/
│   ├─ bvp6c.m
│   └─ (other bvp6c helpers)
├─ InitialShapes/
│   ├─ SIM_Node_50_72_0_1_1_+00_+00.mat
│   ├─ SIM_Node_60_72_0_1_1_+00_+00.mat
│   ├─ SIM_Node_75_72_0_1_1_+00_+00.mat
│   └─ SIM_Node_90_72_0_1_1_+00_+00.mat
├─ SimResults/
│   ├─ solutions
│   │   ├─ 0d7c...5a37.mat
│   │   ├─ 0e17...2304.mat
│   │   └─ (other hash solutions)
│   ├─ cache.mat
│   ├─ catalog.mat
│   └─ OPfile.txt
├─ src/
│   ├─ tools (misc scripts)
│   ├─ utils (helper functions)
│   └─ sim_quad_tree.m
├─ script_driver.mlx
└─ tool_driver.mlx
```