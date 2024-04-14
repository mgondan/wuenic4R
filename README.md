# wuenic4R

Assuming that R.exe is on the PATH,

For a single country
````
cd wuenic4R
estimate.bat bgd
sha1sum.bat bgd
````
The output file is found in the folder out, bgd.txt.

For all countries
````
cd wuenic4R
all.bat
sha1all.bat > v4R.txt
````

Compare SHA to Version 4:
````
diff v4R.txt ..\wuenic4R\v4R.txt
````
