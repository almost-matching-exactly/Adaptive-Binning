Simulation Types:
You can ignore the difference between sim_ and conf_sim_
For all of these, there are 2 covariates.
The simulation type (constraints below) determines the treatment effect. 
All of these below are implicitly multiplied by some beta_tilde
x1 and x2 are both ~Uniform(0, 5)
Constant: everyone 
Box 1D: of the form a < x1 < b
Open box 1D: of the form x1 < b or x1 > a
Box 2D: of the form a < x1 < b and a' < x2 < b'
Open box 2D: of the form a < x1 and a' < x2 (or similar) 
Union box 1D: of the form a < x1 < b or a < x1 < b'
Union box 2D: of the form (a < x1 < b and c < x2 < d) or (a' < x1 < b' and c' < x2 < d')
Linear 1D: of the form a * x1
Linear 2D: of the form a * x1 + b * x2
Quadratic 1D: of the form a * (x1 - c) ^ 2 
Quadratic 2D: of the form a * (x1 - c) ^ 2 + b * (x2 - d) ^ 2
