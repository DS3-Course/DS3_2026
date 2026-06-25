# Working with RStudio projects

- RStudio Projects help organize a self-contained set of directories, scripts, and data files. 
- For more information, see [here](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects) and [here](https://r4ds.hadley.nz/workflow-scripts.html#rstudio-projects) 

## Creating Projects

Create a new project in RStudio via:  
`File > New Project...`

### Starting a new project will:
- Create a `<PROJECT>.Rproj` file (stores project-specific settings). 
- Create a hidden `.Rproj.user` directory (stores temporary files). 
- Load the project in RStudio. 

>[!IMPORTANT]
>The `<PROJECT>.Rproj` file is a shortcut to open the project directly. 

### Opening an existing project in RStudio will:
- Start a new R session. 
- Set the project directory as the current working directory.
- Restore previously open files in the editor.

>[!NOTE]
>You can work on multiple projects at once - each project runs in its own RStudio session.

---

# Organizing your RStudio Project

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

- Only `/data` and R scripts are required - everything else can be recreated. 
- `/data` should be treated as read-only. 
- Analysis outputs go to `/results` or `/plots` (include version info). 
- R workspace and large RDS files are stored in `/rdata`. 
- Additional directories can added as needed (e.g. `/Archive`). 
- Use a standardized workflow where possible in a main R script (e.g. `analysis_v1.0.R`). 
- Keep your main workflow clean by moving reusable code (e.g. functions, common settings etc) to `helper_functions.R` script. 

---

# Using this Project template

## Clone the entire course repository:
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
 
