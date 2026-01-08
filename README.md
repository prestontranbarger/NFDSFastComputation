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

See the Tranbarger paper (at https://arxiv.org/abs/2512.17139) for relevant context.

This section will explain everything you need to know about how to compute the $h_{\gamma,\chi_1,\chi_2,k}$ functions detailed in the paper, in addition to how to graph the $\widehat{S}_{\chi_1,\chi_2,k}$ function to demonstrate quantum modularity.

### 3.1 - Computing $h_{\gamma,\chi_1,\chi_2,k}$

For fixed $k$ and primitive $\chi_1$ and $\chi_2$ with $\chi_1\chi_2(-1)=(-1)^k$ it is simple to compute $h_{\gamma,\chi_1,\chi_2,k}$ for $\gamma\in\Gamma_0(q_1q_2)$.

1. After cloning this repository and creating a new python file within the folder be sure to go through all of the steps in the previous section to make sure your development enviroment is properly set up and that everything follows the expected behaviour.
2. To compute $h_{\gamma,\chi_1,\chi_2,k}$ be sure to include ```from hFunctionInterpolation import *``` in the header of your python file where you would like to compute this function.
3. Consult the example below.

### 3.2 Example $h_{\gamma,\chi_1,\chi_2,k}$ Computation

Given the unique pair of Dirichlet characters $\chi_1$ and $\chi_2$ with moduli $5$ and $5$ respectively such that $\chi_1(2)=\zeta_4^2$ and $\chi_2(2)=\zeta_4^2$ suppose we wish to compute
```math
h_{\gamma,\chi_1,\chi_2,4}\qquad\text{where}\qquad\gamma=\begin{pmatrix}
   51 & 104 \\ 25 & 51
\end{pmatrix}.
```
We can use the following code to accomplish this,
```
from hFunctionInterpolation import *
gamma = matrix(ZZ, [[51, 104], [25, 51]])
dChar1 = allDCharacters(5)[2]
dChar2 = allDCharacters(5)[2]
print(dCharString(dChar1), dCharString(dChar2))
print(hFunctionInterpolate(gamma, dChar1, dChar2, 4))
```

It produces the output,
```
5c5;2-2 5c5;2-2
100%|██████████| 3/3 [00:35<00:00, 11.70s/it]
-24/5*x^2 - 96/5*x - 96/5
```

Which tells us that
```math
h_{\gamma,\chi_1,\chi_2,4}=-\dfrac{24}{5}\,x^2-\dfrac{96}{5}\,x-\dfrac{96}{5}\qquad\text{where}\qquad\gamma=\begin{pmatrix}
   51 & 104 \\ 25 & 51
\end{pmatrix}.
```

### 3.3 - Section 6.2, Precomputed $h_{\gamma,\chi_1,\chi_2,k}$ for the image $\widetilde{S}_{\chi_1,\chi_2,k}(\Gamma_1(q_1q_2))$

Section 6.2 of this paper makes reference to precomputed $h_{\gamma,\chi_1,\chi_2,k}$ which are able to characterize the image $\widetilde{S}_{\chi_1,\chi_2,k}(\Gamma_1(q_1q_2))$ via Theorem 1.17 (for some pairs of characters and weights).

These polynomials can be found in the "hFsOut" folder within the repository, and they are calculated for all quadratic primitive $\chi_1$ and $\chi_2$ with $\chi_1\chi_2(-1)=(-1)^k$. These triples are first indexed by $\chi_1$ and $\chi_2$ by subfolders of the form
```(q_1q_2);(dCharString(dChar1);(dCharString(dChar2)```
Then within these subfolders are further subfolders indexing the weight as follows
```hF-k(weight),(dCharString(dChar1));(dCharString(dChar2))```

Opening one of these files as an example, lets look at "hF-k4,7c7;3-3;3c3;2-1.txt":
```1_1_0_1:0
-20_1_-21_1:144/7*x^2 - 96/49*x
43_-4_441_-41:9072*x^2 - 11808/7*x + 3840/49
85_-9_189_-20:-27270/7*x^2 + 40596/49*x - 2160/49
-146_17_-189_22:14352/7*x^2 - 23124/49*x + 1326/49
22_-3_147_-20:2142*x^2 - 4104/7*x + 1962/49
-146_23_-273_43:-68610/7*x^2 + 150876/49*x - 11856/49
169_-30_231_-41:59340/7*x^2 - 147624/49*x + 13116/49
-188_35_-231_43:14652/7*x^2 - 38184/49*x + 3552/49
85_-16_441_-83:9828*x^2 - 25920/7*x + 17088/49
-62_13_-105_22:-14502/7*x^2 + 42528/49*x - 4458/49
-167_38_-189_43:22710/7*x^2 - 72600/49*x + 8286/49
-209_49_-273_64:42354/7*x^2 - 139008/49*x + 16290/49
106_-25_441_-104:17766*x^2 - 58680/7*x + 48450/49
148_-39_315_-83:-78984/7*x^2 + 291324/49*x - 38382/49
190_-51_231_-62:1464/7*x^2 - 5328/49*x + 684/49
43_-12_147_-41:5292*x^2 - 20736/7*x + 20304/49
22_-7_63_-20:-20166/7*x^2 + 87480/49*x - 13578/49
-125_46_-231_85:-43602/7*x^2 + 224148/49*x - 41136/49
169_-64_441_-167:-49770*x^2 + 263796/7*x - 349536/49
64_-25_105_-41:-17190/7*x^2 + 93972/49*x - 18336/49
-83_34_-105_43:6078/7*x^2 - 34836/49*x + 1020/7
64_-27_147_-62:-3402*x^2 + 20088/7*x - 29646/49
-146_67_-231_106:-130962/7*x^2 + 841800/49*x - 27606/7
-230_107_-273_127:71670/7*x^2 - 466356/49*x + 108384/49
-188_89_-357_169:-104040/7*x^2 + 689520/49*x - 163200/49
211_-100_441_-209:-22680*x^2 + 150480/7*x - 249600/49
232_-123_315_-167:90486/7*x^2 - 671568/49*x + 178014/49
169_-93_189_-104:8640/7*x^2 - 66756/49*x + 18426/49
85_-48_147_-83:-3402*x^2 + 26892/7*x - 53136/49
43_-28_63_-41:18726/7*x^2 - 168468/49*x + 54156/49
106_-75_147_-104:3780*x^2 - 37368/7*x + 92340/49
127_-108_147_-125:1638*x^2 - 19476/7*x + 57888/49
```

On the left are matrices in $\Gamma_1(q_1q_2)$ (which form a finite generating set) written in the form a_b_c_d corresponding to the matrix $\begin{pmatrix}a & b \\ c & d\end{pmatrix}$. Then on the right, after the colon, is the corresponding $h_{\gamma,\chi_1,\chi_2,k}$ for each $\gamma$ in this finite generating set.

### 3.4 - Graphing $\widehat{S}_{\chi_1,\chi_2,k}$

To be added soon.