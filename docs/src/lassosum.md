# lassosum

Mak et al.(2018), Polygenic scores via penalized regression on summary
statistics


Consider linear regression (GWAS) problem

```math
\mathbf{y} = \mathbf{X\beta} + \mathbf{\epsilon}
```

where $\mathbf{X}$ is the n x p genotype matrix (normalized to sd(cols) = 1)
and $\mathbf{y}$ is n x 1 phenotype matrix (normalized to sd(cols) = 1)



Using regularization, we arrive at the LASSO with objective function

```math
f(\mathbf{\beta}) = (\mathbf{y} - \mathbf{X\beta})^T (\mathbf{y} - \mathbf{X\beta}) + 2 \lambda ||\beta||_1
```

```math
f(\mathbf{\beta}) = \mathbf{y}^T \mathbf{y} - 2 \mathbf{\beta}^T \mathbf{X}^T \mathbf{y} + \mathbf{\beta}^T \mathbf{X}^T \mathbf{X} \mathbf{\beta} + 2 \lambda ||\beta||_1
``` 

 
Observing that 

```math
\mathbf{r} = \mathbf{X}^T \mathbf{y}
```

is the genotype - phenotype correlation (obtained from p-value in
summary statistics) and

```math
\mathbf{R} = \mathbf{X}^T \mathbf{X}
```

is the LD matrix, 

the objective can be rewritten as


```math
f(\mathbf{\beta}) = \mathbf{y}^T \mathbf{y} - 2 \mathbf{\beta}^T \mathbf{r} + \mathbf{\beta}^T \mathbf{R} \mathbf{\beta} + 2 \lambda ||\beta||_1
``` 


The LD matrix $\mathbf{R}$ could also be derived from reference
data. In this case, $\mathbf{r}$ and $\mathbf{R}$ are not derived from
the same genotype matrix $\mathbf{X}$ anymore, implying that the
objective function is no longer a penalized least-squares (or LASSO)
problem.

It can be turned into a LASSO problem when regularization

```math
\mathbf{R_s} = (1-s)\mathbf{R} + s\mathbf{I}, 0 \le s \le 1
```

is applied (see simple proof in paper). The shrinkage parameter $s$
determines how much of the linkage disequilibrium between different
variants is taken into account.

This renders the objective function an elastic net,

```math
f(\mathbf{\beta}) = \mathbf{y}^T \mathbf{y} - 2 \mathbf{\beta}^T \mathbf{r} + (1-s)\mathbf{\beta}^T \mathbf{R} \mathbf{\beta} + s \mathbf{\beta}^T \mathbf{\beta} + 2 \lambda ||\beta||_1
``` ,

which can be solved using coordinate descent.

In each iteration of coordiante descent, the contribution of each
single coordinate is optimized on its own. For a single coordinate
$\beta_l$, the contributions to $f(\beta)$ are

```math
\beta_l^2 (s + (1-s)\mathbf{R}_{ll}) - 2 \beta_l (r_l - (1-s)\mathbf{R}_{l}\beta + (1-s)\mathbf{R}_{ll}\beta_l) + 2\lambda |\beta_l|
```


NOTE: Is not $s + (1-s)\mathbf{R}_{ll} = 1$? Only if $R_{ll}$ = 1, but for some markers the variance is 0, hence $R_{ll}$ = 0

Defining 

```math
u_l = r_l - (1-s)R_{.l} \beta + (1-s) R_{ll} \beta_l
```

the minimum contribution is obtained at

```math
\beta_l = \mathmr{sgn}(u_l) (|u_l| - \lambda) / (s + (1-s)\mathbf{R}_{ll}) , |u_l| - \lambda > 0
```

```math
\beta_l = 0, \mathrm{otherwise}
```


NOTE: Is this not the same as is done in LDpred? I.e. we look at the
contribution of all SNPs except $l$ to $\mathbf{r}_l$ due to LD, and
then apply some shrinkage to $\beta_l$.

NOTE: Implementation where thresholded LD matrix is used suffers from
divergence issues. Implementation using genotype matrix directly
converges fine.