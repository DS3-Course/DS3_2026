# HTP DESeq2 analysis

**Temporary fix for biobroom tbl_df error:**
assignInNamespace("tbl_df", tibble::as_tibble, ns = "dplyr")
add this to your main R script.

1. **Ensure you have the course materials** 
      1. **Clone** this repository (first time only)
         ```bash
         git clone https://github.com/DS3-Course/DS3_2026.git
         ```
      2. **Pull** updates each day to ensure updated content:
         ```bash
         git pull
         ```
         *Alternative:* direct download from GitHub

---

2. **Open the project in RStudio**  
      - **Option A:** Navigate to the folder and double-click the `.Rproj` file  
      - **Option B:** In RStudio, go to `File > Open Project...` and select the `.Rproj` file  

> 💡 Always open the `.Rproj` file first (not individual scripts)
  
> ⚠️ This project uses `renv` for package management (install with `install.packages("renv")` if needed)

---

3. **Run the analysis**

   - Open the R script(s)  
   - Work through the analysis steps in the script  