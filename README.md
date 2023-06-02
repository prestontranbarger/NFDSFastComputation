# NFDSFastComputation

## 1 - Windows Installation & Setup

### 1.1 - Cloning the Repository
1. Install git from this URL ```https://git-scm.com/download/win```; most machines should be 64-bit windows. Use the default selections for the installer.
2. Using file explorer, create a new directory in a convenient location on your machine; we will clone this repository into this folder in a later step.
3. Open Git Bash and navigate to the new directory you just created using the ```cd``` and ```ls``` commands.
4. Now clone this repository using the command ```git clone https://github.com/prestontranbarger/NFDSFastComputation```. If you wish to copy-paste this command into Git Bash you must right-click and select paste.

### 1.2 - Installing Sage and Setting up the IDE
1. Download the SageMath 9.1 Windows installer from this URL ```https://github.com/sagemath/sage-windows/releases```.
2. Once this finishes, add the following to your Path enviroment variables (here the ellipses denote the location of your SageMath 9.1 installation, more information can be found here ```https://ask.sagemath.org/question/39742/make-pycharm-recognise-the-sage-python-interpreter/```):
   - ```...\SageMath 9.1\runtime\opt\sagemath-9.1\local\lib```
   - ```...\SageMath 9.1\runtime\opt\sagemath-9.1\local\bin```
   - ```...\SageMath 9.1\runtime\bin```
   - ```...\SageMath 9.1\runtime\lib\libpack```
3. Now create the following new enviroment variables:
   - DOT_SAGE: ```/dot_sage```
   - SAGE_LOCAL: ```/opt/sagemath-9.1/local```
   - SAGE_ROOT: ```/opt/sagemath-9.1```
4. Download the PyCharm community version installer from this URL ```https://www.jetbrains.com/pycharm/download/#section=windows```. Use the default selections for the installer.
5. Once PyCharm is downloaded, create a new project under ```File >New Project```. All of the setup instructions for this project are in the following steps (do not immediately click ```Create```).
6. For the ```Location``` parameter of our project, first navigate to the directory you created in section 1.1 step 2; then select the subdirectory with title ```NFDSFastComputation```.
7. For the ```Python Interpreter``` parameter of our project select ```Previously Configured Interpreter```; then select the three dots to the right of the ```Interpreter``` dropdown.
8. Select ```System Interpreter```; then select the three dots to the right of the ```Interpreter``` dropdown.
9. Navigate to the location of your SageMath 9.1 installation (you may need to click a button to allow hidden folders depending on your installation location).
10. Now naviage to the following subdirectory ```runtime\opt\sagemath-9.1\local\bin``` and select ```python3.7.exe```. Select ```OK``` to finalize the location of your interpreter.
11. Select ```OK``` to add the SageMath 9.1 interpreter; then select ```Create``` and ```Create from Existing Sources``` to open the project in the IDE.
