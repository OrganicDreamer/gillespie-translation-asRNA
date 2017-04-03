# To use with PyCharm/IDE

1) Download and install <a href='https://www.jetbrains.com/pycharm/download/#section=windows'>PyCharm</a>
2) <a href='https://www.continuum.io/downloads'>Download</a> and <a href='https://docs.continuum.io/anaconda/install'>install</a> Anaconda. 
3) <a href='https://docs.continuum.io/anaconda/ide_integration#pycharm'>Set up PyCharm</a> to use Anaconda.
4) Download and unzip this folder
5) Open .py files in this folder in Pycharm. Ensure .py scripts are opened with this folder as the directory.
6) Use CTRL-A to select all and ALT-SHIFT-E to run script of choice in console

# To run headless/on server cluster (specifically for Windows)

1) Download appropriate Windows installer for <a href ='http://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html'>PuTTY</a> 
2) Download <a href='https://winscp.net/eng/download.php>WinSCP</a> installer
3) Connect and log in to server using both
4) Download and unzip this folder
5) Using WinSCP,copy/move the files in the folder onto home directory on cluster
6) Using PuTTY, run bash script from home directory on the command line (qsub runme.sh)
