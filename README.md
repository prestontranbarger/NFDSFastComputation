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

## 2 - Fast Computation of Generalized Dedekind Sums
See the Tranbarger, Wang paper (at https://www.worldscientific.com/doi/10.1142/S179304212450060X) for relevant context.

This section will explain everything you need to know about how to leverage the algorithm detailed in the above paper to compute weight $2$ generalized Dedekind sums in an efficient manner. For the higher weight case, continue to the next section.

### 2.1 - Using the Repository
1. Create a new python file in the ```NFDSFastComputation``` directory in your PyCharm enviroment.
2. To use the NFDS methods in this python file, include ```from  NFDS import *``` at the top of your code. The main method you will be using is ```newFormDedekindSumFast(dChar1, dChar2, gamma, chpr)```. Below we have example code for how to use this method; however, you will also be using methods such as ```allDCharacters(n)```, ```modulus(dChar)```, ```isEven(dChar)```, ```isPrimitive(dChar)```, ```chprPathFinder(dChar1, dChar2)```, and ```readAllChpr(path)```. Please email me with questions about how to use these methods if there is confusion.
3. When running your code for the first time you may need to install some packages, use the ```Python Packages``` tab near the bottom of your IDE to install these packages. For example, if you recieve the error message ```ModuleNotFoundError: No module named 'tqdm'``` then use the package manager to install the ```tqdm``` library.

### 2.2 - Example
Given the unique pair of Dirichlet characters $\chi_1$ and $\chi_2$ with moduli $5$ and $7$ respectively such that $\chi_1(2)=\zeta_4^3$ and $\chi_2(3)=\zeta_6^5$ suppose we wish to compute
```math
S_{\chi_1,\chi_2}\begin{pmatrix}37&18\\35&17\end{pmatrix}.
```
We can use the following code to accomplish this,
```
from NFDS import *

dChar1 = allDCharacters(5)[3]
print(dChar1(2))
dChar2 = allDCharacters(7)[5]
print(dChar2(3))

print(isEven(dChar1) == isEven(dChar2), isPrimitive(dChar1), isPrimitive(dChar2))
print(dCharString(dChar1), dCharString(dChar2))

mtrx = Matrix([[37, 18], [35, 17]])

print(newFormDedekindSumFast(dChar1, dChar2, mtrx, readAllChpr(chprPathFinder(dChar1, dChar2))))
```

It produces the output,
```
-zeta4
-zeta6 + 1
True True True
5c5;2-3 7c7;3-5
1.4928203230278922 + 0.9856406460553006*I
100%|██████████| 1/1 [00:00<00:00, 539.67it/s]
```

Which tells us that our Dedekind sum equals the resulting complex number.

### 2.3 - Notes
- When handling the raw precomputation data for pairs of characters, it is best to use the ```characterPairsLM``` data as this is the most expansive data set in the repository.

## 3 - Higher Weight Generalized Dedekind Sums
See the Tranbarger paper (at https://coming.soon) for relevant context.

TO BE ADDED WITH ARXIV RELEASE