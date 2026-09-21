# An Adjustable Vector Topology Analysis Method for the Magnetic Gradient Tensor Field
## Basic Usage and Installation

No installation or compilation of the source code is required. The programs in this repository were developed and tested using **MATLAB R2020b**.

The repository is organized by figure. Each `fig` folder is a **self-contained unit corresponding exclusively to the same-numbered figure in the paper**. For example, the files in the `fig4` folder are provided specifically for reproducing Fig. 4 and should be used within the `fig4` folder. Files from different `fig` folders should not be moved between folders or mixed together.

Each `fig` folder already contains:

* the MATLAB `.m` files required for the corresponding calculation and visualization; and
* the input data required by these programs, provided as `.mat` files.

The commands for loading the required input data have already been included in the corresponding MATLAB programs. Therefore, users do not need to manually import the data, specify the input data files, or modify the data-loading commands when the original folder structure is retained.

To reproduce a result:

1. Install **MATLAB R2020b** and the required MATLAB toolbox listed in the **Dependencies and Computational Requirements** section.
2. Clone the repository or obtain the complete `fig` folder corresponding to the figure to be reproduced.
3. Keep the contents of that `fig` folder unchanged. Do not move files from one `fig` folder to another.
4. Open MATLAB and set the corresponding `fig` folder as the current working directory.
5. Run the `.m` file whose filename ends with `_main` or `_draw`.
6. The program will automatically use the input data provided in the same folder and perform the corresponding calculations and/or generate the visualization.

For example, to reproduce Fig. 4, use the programs and data contained in the `fig4` folder and run the corresponding `_main.m` or `_draw.m` file from that folder. The files in `fig4` are intended for Fig. 4 and should not be used as input files for programs in other `fig` folders.

## Dependencies and Computational Requirements

### Dependencies

The programs were developed and tested using the following software:

* **MATLAB R2020b**
* **Statistics and Machine Learning Toolbox**, required for the `dbscan` function used in the programs.

All input data required by the programs are already provided in their corresponding `fig` folders, and the commands for loading these data are included in the programs. Therefore, no additional data files or data-processing software are required.

### Computational Requirements

No specialized computing hardware or high-performance computing environment is required. The programs can be executed on a standard desktop or laptop computer capable of running MATLAB R2020b and the required toolbox.

No GPU computation is required.


## License

This project is licensed under the Apache License 2.0. See the `LICENSE` file for details.
