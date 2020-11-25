# UNC Bioinformatics Course

Welcome! Please refer to the Course Outline and Course Schedule for more information on delivery and assessment.

# Software requirements

Please ensure the following programs are installed on your desktop **prior to commencing the class**. 

* R https://www.r-project.org/
* R-Studio https://rstudio.com/
* Cisco AnyConnect VPN https://ecu.teamdynamix.com/TDClient/1409/Portal/Requests/ServiceDet?ID=11945 & https://itcs.ecu.edu/2020/03/16/working-remotely/

Some method for accessing the ECU computers via unix commandline is required. 

## Windows Users

Windows 10 offers a linux subsystem that enables you to use linux scripts. https://docs.microsoft.com/en-us/windows/wsl/install-win10. 
RECOMMENDED: I use the Ubuntu subsystem (Step 6), and many bio servers are based off Ubuntu so this is the best option.

If you are not running Windows 10, you can use PuTTY to access the ECU biology server (https://www.putty.org/). Please contact me if you are struggling to set up this option.

## Mac Users

You already have access to unix scripting through your OSX terminal. You can search for it using "Terminal" or find it in the ```Utilities``` folder in ```Applications```


## Required R packages 

The following packages need to be installed on R on your desktop. 
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("rhdf5")