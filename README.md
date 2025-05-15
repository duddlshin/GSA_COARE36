# SensitivityAnalysis_COARE

Repository for project in MIT class 12.870 (Air-Sea Interactions) in which a global sensitivity analysis is conducted for empirical coefficients in the COARE 3.6 algorithm (https://github.com/NOAA-PSL/COARE-algorithm.git).

### Target coefficients for sensitivity analysis:
- $m$ and $b$, which formulate the Charnock coefficient $\alpha = m U_{10N} + b$ where $U_{10N}$ is the wind speed at reference height.

- roughness Reynolds number $\gamma$.

### Variance-based method for sensitivity analysis:
- Two metrics are used in a variance-based analysis

    1. First-order sensitivity indices: how much variance in the model output $Y$ is explained by the coefficient $X_i$ alone, that is, discluding its interacting effects with other coefficients.

        $S_i = \frac{Var_{X_i} \( E_{X_{\sim i}}[Y \mid X_i] \) }{Var(Y)}$
 
    2. Total sensitivity indices: how much variance in the model output is explained by the overall effect of coefficient $X_i$, inclduing its interacting effects.
 
        $S_{T_i} = \frac{ E_{X_{\sim i}} \(Var_{X_{i}} [Y \mid X_{\sim i}] \) }{Var(Y)}$

- Both metrics are estimated using the method in Jansen, 1999 (https://www.sciencedirect.com/science/article/abs/pii/S0010465598001544). The steps are as follow:

    1. Two matrices of size ($N \times m$) are initialized where $N$ is the number of samples and $m$ is the number of coefficients. We will call these matrix $A$ and $B$.
    2. These are propagated forward into the model to get response matrices $f(A)$ and $f(B)$ that contain the model outputs given the columns of coefficients.
    3. A single column, i.e., single coefficient, in matrix $A$ is switched with the corresponding column in matrix $B$. Effectively, a single coefficient is varied, while others are kept constant. This is done for each parameter and run forward into the model. Therefore, the method requires $N(2+m)$ evaluations of the algorithm.
    4. These two model output matrices are used to estimate the components in the first-order and total sensitivity indices:

        $Var_{X_i} \left( E_{X_{\sim i}}[Y \mid X_i] \right) \approx Var(Y) - \frac{1}{2N} \sum_{j=1}^{N} \left( f(B)_j - f(A_B^{(i)})_j \right)^2 $
 
        $E_{X_{\sim i}} \(Var_{X_{i}} [Y \mid X_{\sim i}] \) \approx \frac{1}{2N} \sum_{j=1}^{N} \left( f(A)_j - f(A_B^{(i)})_j   \right)^2 $


### How to use the code
The code for COARE 3.6 algorithm is modified to accept random values for the target coefficients. 
Step through $\textit{main.m}$ code to run the sensitivity analysis.
The only values that need to be configured are the sample size $N$ used to compute the indices (found in line 16 of $\textit{main.m}$) and the model output (found in line 17 of $\textit{main.m}$) for which the sensivity to coefficients is being assessed.

All questions, comments, and feedback are welcome and may be addressed to: youngin@mit.edu.
