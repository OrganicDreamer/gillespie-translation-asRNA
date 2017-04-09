# Set up PyCharm/IDE

1) Download and install <a href='https://www.jetbrains.com/pycharm/download/#section=windows'>PyCharm</a> or your IDE of choice.
2) <a href='https://www.continuum.io/downloads'>Download</a> and <a href='https://docs.continuum.io/anaconda/install'>install</a> Anaconda. 
3) <a href='https://docs.continuum.io/anaconda/ide_integration#pycharm'>Set up PyCharm</a> to use Anaconda.

# To use with PyCharm/IDE
4) Download and unzip this folder
5) Open .py files in this folder in Pycharm/IDE.
6) Use CTRL-A to select all and ALT-SHIFT-E to run script of choice in console

# Set up to run headless/on server cluster (specifically for Windows)

1) Download <a href='https://winscp.net/eng/download.php>WinSCP</a> installer
2) Connect and log in to server using both
3) Download and unzip this folder
4) Using WinSCP,copy/move the files in the folder onto home directory on cluster

# To run headless/on server cluster (specifically for Windows)

1) Download appropriate Windows installer for <a href ='http://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html'>PuTTY</a> 
2) Open PuTTY and log in to the server.
3) Copy/modify .py files to model desired simulation parameters.
4) Copy/modify setup_run_py_script.sh so it targets the desired .py simulation(s) you want to run.
5) Modify RUN_THIS.sh so it will add different setup_run_py_script.sh to run on the cluster simulataneously, if desired. 
6) Add bash script to queue from home directory on the command line (qsub RUN_THIS.sh). Command line: qstat to monitor jobs on the cluster.
