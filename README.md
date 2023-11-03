[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=mathworks/open-gravitational-wave-data-matlab&file=main/open-gravitational-wave-data-matlab/blob/main/OpenPhysicsTutorial.mlx)

# Analyse Open Gravitational Wave Data in MATLAB®

A MATLAB Live Script with accompanying 
- Jupyter® Notebook, 
- m file and 
- reproducible code capsule on Code Ocean®

to access and analyze Gravitational Wave data sets from the **Gravitational Wave Open Science Center (GWOSC)** database

## Get started

Use this tutorial to get started with freely available gravitational wave data at [GOWSC](https://gwosc.org/) directly from MATLAB.
- No downloads, no installations
- **Open directly in MATLAB Online** by clicking this ![badge](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)
- Step-by-step tutorial shows how to
    - Re-use available gravitational wave data. **Access a list of openly available events and sessions** from LIGO, Virgo and Kagra
    - **Query and inspect the metadata** associated with these projects using commands directly from MATLAB (RESTful API)
    - Avoid downloads. **Access specific data** from within the database directly and **avoid time-consuming downloads** of large data
    - **Read in HDf5 data** from within MATLAB corresponding to specific sampling rates and detectors
    - **Analyze gravitational wave data** using standard signal processing in time and frequency domains **re-using [open community](https://de.mathworks.com/matlabcentral/fileexchange/108859-gravitationalwavedataexplorer?tab=example&focused=) code from File Exchange**
    - Let others run your code and reproduce your results quickly. **Pubish the results on GitHub** and **make them accessible** using Open With MATLAB Online
    - Allow people to cite you! **Generate a DOI** for your code by linking your GitHub repository to one of several DOI-generating sites.
- **Live Script** contains **easy-to-use menus** for user to click and select different datasets
- Available on [File Exchange](mathworks.com/matlabcentral/fileexchange/) for directly installing onto your MATLAB path with one click using the [Add-Ons button](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html)
- Accompanying **Jupyter notebook** (.ipynb) for use in a Jupyter environment. More information on MATLAB kernel [here](mathworks.com/products/reference-architectures/jupyter.html)
- Accompanying **Code Ocean reproducible capsule** (.m) for one-click reproducibility of the code by anyone, including reviewers.

## About the Gravitational Wave Open Science Center 
The Gravitational Wave Open Science Center (GWOSC) is a public repository of gravitational wave events and experimental sessions from LIGO, Virgo and Kagra for the community.
It can be accessed at [https://gwosc.org/](https://gwosc.org/). [Here](https://gwosc.org/eventapi/html/GWTC/) is a list of available GW events from GWOSC.

**For advanced users** A detailed guide to the GWOSC API can be found [here](https://gwosc.org/apidocs/). To access the REST API use the MATLAB [webread](mathworks.com/help/matlab/ref/webread.html) function

### Required Products
This tutorial uses the following products\
- MATLAB
- Signal Processing Toolbox

This code has been developed and tested using MATLAB 2023a

**Note**
This tutorial works best when delivered by a tutor. It is important to highlight best practices when working with Open Data, publishing Open Code or making research output reproducible
