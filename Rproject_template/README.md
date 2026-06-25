# Working with RStudio projects

- RStudio projects make it easier to create and manage a self-contained set of directories, scripts, and data files. 
- For more information, see [here](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects) 

## Creating Projects
To create a new Project in  RStudio, go to  `File > New Project...`  

Starting a new project:
- Creates a `<PROJECT>.Rproj` file in the project directory. This file stores project-specific options. 
- Creates a hidden directory `.Rproj.user` that stores project-specific temporary files. 
- Loads the project in RStudio. 

>[!IMPORTANT]
>The `<PROJECT>.Rproj` file can be used as a shortcut to open the project directly from your filesystem. 

Opening a project in RStudio (among other things):
- Starts a new R session. 
- Sets the project directory as the current working directory.
- Previously open source documents are restored in editor tabs.

>[!NOTE]
>You can work on multiple projects at once - each project runs in its own instance of RStudio.


## Organizing your RStudio Project

```
Project_directory
├ /data
├ /results
├ /plots
├ /rdata
├ helper_functions.R
├ analysis_script.R
├ project.Rproj
```

- Only `/data` and R scripts are required - everything else can be recreated (incl. earlier versions). 
- `/data` should be treated as read-only. 
- Analysis outputs go to `/results` or `/plots` (with version info). 
- R workspace and large RDS are stored in `/rdata`. 
- Additional directories can added as needed, eg `/Archive`. 
- Use a standardized workflow where possible in a main `<analysis_vN.N>.R` script. 
- Keep your main workflow clean by moving functions, common settings etc to `helper_functions.R` script. 

## Using this Project template

### Clone the entire course repository:
1. **Clone** this repository (first time only)
   ```bash
   git clone https://github.com/DS3-Course/DS3_2026.git
   ```
2. **Pull** updates each day to ensure updated content:
   ```bash
   git pull
   ```
   *Alternative:* direct download from GitHub
3. **Open** the appropriate .Rproj file for each lesson in RStudio
4. Open and edit R scripts
 
