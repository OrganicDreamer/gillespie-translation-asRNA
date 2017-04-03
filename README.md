# To use with PyCharm/IDE

1) Download and unzip this folder
2) Open .py files in this directory in Pycharm and select Python Interpreter (Anaconda recommended)
3) Use CTRL-A to select all and ALT-SHIFT-E to run script of choice in console

# To run headless/on server cluster (specifically for Windows)

1) Download appropriate Windows installer for <a href ='http://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html'>PuTTY</a> 
2) Download <a href='https://winscp.net/eng/download.php>WinSCP</a> installer
3) Connect and log in to server using both
4) Download and unzip this folder
5) Using WinSCP,copy/move the files in the folder onto home directory on cluster
6) Run bash script (qsub thefile.sh)
