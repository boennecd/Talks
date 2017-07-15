dynamichazard Dynamic Hazard Models using State Space Models
========================================================
author: Benjamin Christoffersen
date: 28/07/2017
autosize: true
css: custom.css

<script>
var h1s = document.getElementsByTagName("h2");
h1s[0].outerHTML = "<h1 style='margin-bottom: 0;'>dynamichazard</h1><h2 style='margin-top: 0;'>Dynamic Hazard Models <br /> using State Space Models</h2>";
//(function($){ 
//  $( "h1" ).html( "Dynamic Hazard Models <br /> using State Space Models: dynamichazard" ); 
// })(jQuery); 
</script> 

<!-- Two line header -->
<script>
var h1s = document.getElementsByTagName("h1");
h1s[0].outerHTML = "<h1 style='margin-bottom: 0;'>dynamichazard</h1><h2 style='margin-top: 0;'>Dynamic Hazard Models <br /> using State Space Models</h2>";
//(function($){ 
//  $( "h1" ).html( "Dynamic Hazard Models <br /> using State Space Models: dynamichazard" ); 
// })(jQuery); 
</script>  

<!-- See this link for macros with mathjax: http://docs.mathjax.org/en/latest/tex.html#tex-macros -->
<script>
MathJax.Hub.Config({
  TeX: {
    Macros: {
      mat: ["\\mathbf{#1}", 1],
      
      Lparen: ["\\left( #1\\right)", 1],
      Lbrace: ["\\left\\{ #1\\right\\}", 1],
      
      Cond: ["\\left. #1 \\vphantom{#2} \\right\\vert #2", 2],
      
      emNotee: ["#1_{\\left. #2 \\right\\vert #3}", 3],
      
      Prob: ["\\text{P}"],
      propp: ["\\Prob\\Lparen{#1}", 1],
      proppCond: ["\\propp{\\Cond{#1}{#2}}", 2],
      
      E: ["\\text{E}"],
      expecp: ["\\E\\Lparen{#1}", 1],
      expecpCond: ["\\expecp{\\Cond{#1}{#2}}", 2],
      
      deter: ["\\left| #1 \\right|", 1],
      
      argminu: ["\\underset{#1}{\\text{argmin}}\\:", 1]
    }
  }
});
</script>




Presentation
========================================================

* Model
* Hard disk failures
* Estimation methods
* End comments

Model
========================================================

Individual $i = 1,\dots,n$

Failure times $T_1,\dots,T_n\in \mathbb{R}_+$

Timer intervals $t = 1, \dots, d$

Binary outcomes $y_{it} = \left\{ \begin{matrix} 1 & T_i \in (t-1, t] \\ 0 & \text{otherwise} \end{matrix}\right.$

Censoring indicators $D_{i1},\dots,D_{i2},\dots$ <small> (one if censored)</small>

Covariates $\vec{x}_{i1}, \vec{x}_{i2}, \dots$

Model
========================================================

## Observational model

$$\proppCond{Y_{it} = 1}{ \vec{\alpha}_t} =  h(\vec{\alpha}_t^\top \vec{x}_{it})$$

$$h(\eta) = 1 / (1 + \exp(-\eta))$$

## State model

$$\begin{array}{ll}
  \vec{\alpha}_{t + 1} = \mat{F}\vec{\alpha}_t + \mat{R}\vec{\eta}_t \qquad &
  \vec{\eta}_t \sim N(\vec{0}, \mat{Q}) \\
  & \vec{\alpha}_{0} \sim N(\vec{a}_0, \mat{Q}_0)
\end{array} , \qquad t = 1,\dots, d$$



Hard disk failures
========================================================
left: 50%



<img class = "icon" alt="BackBlaze logo" src="dl-figure/backblaze-logo.gif"/>

Storage provider

More than 65000 hard disks Today

3 years of data

84366 observations and 623268 data rows.


***
![vault image](dl-figure/backblaze-b2-07-datacenter-aisle-half.jpg)

Hard disk failures
========================================================


```r
hd_dat[
  1:10, c("serial_number", "model", "tstart", 
          "tstop", "fails", "smart_12")]
```

```
    serial_number        model tstart tstop fails smart_12
505      5XW004AJ ST31500541AS   30.0  40.0     0        0
506      5XW004AJ ST31500541AS   40.0  43.2     0       24
507      5XW004AJ ST31500541AS   43.2  56.9     0       25
508      5XW004Q0 ST31500541AS   40.6  51.0     0        0
509      5XW004Q0 ST31500541AS   51.0  53.7     0       54
510      5XW004Q0 ST31500541AS   53.7  54.1     0       56
511      5XW004Q0 ST31500541AS   54.1  54.4     0       57
512      5XW004Q0 ST31500541AS   54.4  54.5     0       58
513      5XW004Q0 ST31500541AS   54.5  54.7     0       59
514      5XW004Q0 ST31500541AS   54.7  57.2     0       61
```

Estimation methods: Log likelihood
========================================================
$$
L(\mat{Q},\mat{Q}_0, \vec{a}_0)
	= p(\vec{\alpha}_0)\prod_{t = 1}^d \proppCond{\vec{\alpha}_t}{\vec{\alpha}_{t-1}}
		\prod_{i \in R_t} \proppCond{y_{it}}{\vec{\alpha}_t}
$$
 
Risk set $R_t = \Lbrace{i \in \mathbb{Z}_{+}: D_{it} = 0}$


Estimation methods: Log likelihood
========================================================

$$
\begin{aligned}
	\log L \Lparen{\mat{Q},\mat{Q}_0, \vec{a}_0} \propto &
	 \mathcal{L}\Lparen{\mat{Q},\mat{Q}_0, \vec{a}_0}
		 \\
		= & - \frac{1}{2} \Lparen{\vec{\alpha}_0 - \vec{a}_0}^\top \mat{Q}^{-1}_0\Lparen{\vec{\alpha}_0 - \vec{a}_0} \\
	&  - \frac{1}{2} \sum_{t = 1}^d \Lparen{\vec{\alpha}_t - \mat{F}\vec{\alpha}_{t - 1}}^\top\mat{R}^\top\mat{Q}^{-1}\mat{R}\Lparen{\vec{\alpha}_t - \mat{F}\vec{\alpha}_{t - 1}} \\
	&  - \frac{1}{2} \log \deter{\mat{Q}_0} - \frac{1}{2d} \log \deter{\mat{Q}} \\
	&  + \sum_{t = 1}^d \sum_{i \in R_t} l_{it}({\vec{\alpha}_t})
\end{aligned}
$$

Estimation methods: Log likelihood
========================================================

## E-step
Find (smoothed) estimates of $\vec{\alpha}_0,\dots, \vec{\alpha}_d$ and smoothed covariance matrices given $\mat{Q}$, $\mat{Q}_0$ and $\vec{a}_0$.

## M-step
Update $\mat{Q}$ and $\vec{a}_0$.

### References
<small>Due to Fahrmeir (1992) and Fahrmeir (1994). Like Shumway and Stoffer (1982) which if for the standard Kalman filter. </small>

Estimation methods: E-step
========================================================

## Filtering
* Extended Kalman Filter
* Unscented Kalman filter
* Sequential approximation of posterior modes
* Mode estimation

## Smoothing
Rauch-Tung-Striebel algorithm

Estimation methods: Kalman Filter (KF)
========================================================

$$
\emNotee{\vec{a}}{t}{s} = \expecpCond{\vec{\alpha}_t}{\vec{y}_1,\dots,\vec{y}_s},  \qquad
    \emNotee{\mat{V}}{t}{s} = \expecpCond{\mat{V}_t}{\vec{y}_1,\dots,\vec{y}_s}
$$

## Prediction step
Compute $\emNotee{\vec{a}}{t}{t - 1}$ &  $\emNotee{\mat{V}}{t}{t - 1}$ with $\emNotee{\vec{a}}{t - 1}{t - 1}$ and $\emNotee{\mat{V}}{t - 1}{t - 1}$

## Correction step
Compute $\emNotee{\vec{a}}{t}{t}$ &  $\emNotee{\mat{V}}{t}{t}$ with $\emNotee{\vec{a}}{t}{t-1}$ and $\emNotee{\mat{V}}{t}{t - 1}$

Estimation methods: Extended Kalman Filter (EKF)
========================================================

Same prediction step as KF

Make 1. order Taylor expansion around $\emNotee{\vec{a}}{t}{t - 1}$ to get similar correction step

Same as one Fisher scoring step in:

$$
\begin{array}{c}
\emNotee{\vec{a}}{t}{t} = \argminu{\vec{\alpha}}
  - \log \proppCond{\vec{\alpha}}{\emNotee{\vec{a}}{t}{t-1}, \emNotee{\mat{V}}{t}{t-1}}
    -\sum_{i \in R_t} \log\proppCond{y_{it}}{\vec{\alpha}} \\
%
\vec{\alpha} \sim N\Lparen{\emNotee{\vec{a}}{t}{t-1}, \emNotee{\mat{V}}{t}{t-1}}
\end{array}
$$

Essentially L2 penalized generalized linear models.

Estimation methods: mode approximation
========================================================

(Global) mode approximation: Minimize likelihood direclty using Newton Raphson 

EKF with extra iteration: Take more iteration. Use working responses as `glm`

Examples
========================================================




References
========================================================
<small><div class="references">
[1] L. Fahrmeir. "Dynamic Modelling and Penalized Likelihood Estimation for Discrete Time Survival Data". In: _Biometrika_ 81.2 (1994), pp. 317-330.
[1] L. Fahrmeir. "Posterior Mode Estimation by Extended Kalman Filtering for Multivariate Dynamic Generalized Linear Models". In: _Journal of the American Statistical Association_
87.418 (1992), pp. 501-509.
</div></small>
