# Getting started with `R` and `RStudio`

## 1. Installing `R`

a. For Mac OS users, you can download `R` from CRAN [here](https://cran.r-project.org/bin/macosx/).  
   - For R ≥4.1.0 there is a separate build for ARM-based Macs (also known as ‘M1/M2/M3/M4’ and ‘Apple Silicon’).  
   - _For packages that require compilation:_ Xcode Command Line Tools needs to be installed - to install execute the following command in Terminal: `xcode-select --install` or see [here](https://cran.r-project.org/bin/macosx/tools/) and [here](https://mac.install.guide/commandlinetools). If having gfortran issues on Apple Silicon computers, check out the [macrtools R package](https://github.com/coatless-mac/macrtools/), or use the gfortran installer linked at the top of [this page](https://cran.r-project.org/bin/macosx/tools/).   

b. For Windows users, you can download `R` [here](https://cran.r-project.org/bin/windows/base/).  
   - _For packages that require compilation:_ see [here](https://cran.r-project.org/bin/windows/Rtools/).  

c. For Linux users, you can download `R` [here](https://cran.rstudio.com/bin/linux/).  

## 2. Installing `RStudio` IDE  
Download `RStudio` [here](https://posit.co/download/rstudio-desktop/)  
  
For tips on getting started, see:  
   - `RStudio` [User Guide](https://docs.posit.co/ide/user/)  
   - `RStudio` [Get Started](https://docs.posit.co/ide/user/ide/get-started/) section.
   - `RStudio` [cheatsheet](https://rstudio.github.io/cheatsheets/html/rstudio-ide.html)  

## 3. Updating `R` and `RStudio`  

For more information on how `RStudio` and `R` interact, see [here](https://support.posit.co/hc/en-us/articles/200486138-Changing-R-versions-for-the-RStudio-Desktop-IDE).  

a. Mac OS:
   - Simply downloading and installing a new version of `R` should result in RStudio using the latest version.  
   - By default `RStudio` will run against the current version of R.Framework.  
   - List all of the versions of R.Framework on your system and determine which one is considered the current one by executing the following command in Terminal:
		```
		ls -l /Library/Frameworks/R.framework/Versions/
		```
   - To change the current version of R.Framework, you can either:  
     a. Run the installer from CRAN for the R version you want to be current.  
     b. Update the R.framework/Versions/Current directory alias directly using `ln -s`.  
	  c. Override the RSTUDIO_WHICH_R environment variable to the R executable that you want to run against (add to ~/.profile file):
		```
		export RSTUDIO_WHICH_R=/usr/local/bin/R
		```

b. Windows:
   - You may need to uninstall `R` to update it.  
   - On Windows, `RStudio` uses the system's current version of `R` by default.
   - You can override which version of `R` is used via General panel of the `RStudio` Options dialog.
   - Another option for Windows users is to use a package called `installr`

c. Linux:
   - Simply downloading and installing a new version of `R` should result in `RStudio` using the latest version.  
   - `RStudio` uses the version of `R` pointed to by the output of the following command:
		```
		which R
		```
   - To override the RSTUDIO_WHICH_R environment variable to the `R` executable that you want to run against (add to ~/.profile file):
		```
		export RSTUDIO_WHICH_R=/usr/local/bin/R
		```

d. Updating `RStudio`
   - Go to `Help > Check for Updates > Quit and Download...`  
   - Download and install.  


## Alternatives to running R and RStudio locally
a. Data Studio on the [CAVATICA platform](https://cavatica.sbgenomics.com/) (free account required)

b. [Posit Cloud](https://posit.cloud/plans/free) (free account required)


## Additional requirements
To allow cloning of project repositories from within RStudio, you will need to have Git and (optionally) Git LFS installed on your system. 

### Git
For details, see [here](https://git-scm.com/downloads)  

a. Mac OS:  
   - Using [Homebrew](https://brew.sh/): `brew install git`   
   - Using [MacPorts](https://www.macports.org/): `sudo port install git`   

b. Windows:  
   - Download installer [here](https://git-scm.com/download/win)  
   - Using [WinGet](https://learn.microsoft.com/en-us/windows/package-manager/winget/): `winget install --id Git.Git -e --source winget` 

c. Linux:  
   - Instructions [here](https://git-scm.com/download/linux)  
   - Debian/Ubuntu: `apt-get install git`  
   - Fedora/Rocky: `dnf install git `

### Git LFS (large file storage)  
OPTIONAL - only required during DS3 course for mass_cytometry_dataviz exercise.  
For details, see [here](https://git-lfs.com/)  

a. Mac OS:  
   - Using [Homebrew](https://brew.sh/): `brew install git-lfs`   
   - Using [MacPorts](https://www.macports.org/): `port install git-lfs`   

b. Windows:  
   - Git LFS is included in the distribution of Git for Windows  

c. Linux:  
   - Instructions [here](https://github.com/git-lfs/git-lfs/blob/main/INSTALLING.md)  

Once installed, set up Git LFS for your user account by running:
`git lfs install`

