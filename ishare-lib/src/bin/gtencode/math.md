### 1. The Rare Variant Mask (Filtering for Rare Variants)

To define this mathematically, we introduce an **indicator function** (or mask) based on the population allele frequency, $p_{m,k}$. Let $\tau$ be your threshold for "rare" (for example, a Minor Allele Frequency of $\tau = 0.01$).

We define the weight/mask $w_{m,k}$ for allele $k$ at locus $m$ as:
$$w_{m,k} = \begin{cases} 1 & \text{if } p_{m,k} \le \tau \\ 0 & \text{if } p_{m,k} > \tau \end{cases}$$

Now, we define the "rare-filtered" allele count for individual $X$ as:
$$\tilde{x}_{m,k} = x_{m,k} \cdot w_{m,k}$$

If an allele is common, $\tilde{x}_{m,k}$ becomes $0$, entirely erasing its contribution to both individuals' vectors. If it is rare, the original count (1 or 2) is preserved.


### 2. Rare-Variant Cosine Similarity

The magnitude of the vector depends *exclusively* on how many rare alleles that specific individual carries.

**Mathematical Definition:**
$$\text{Cosine}_{rare}(X, Y) = \frac{\sum_{m=1}^{M} \sum_{k=1}^{K_m} \tilde{x}_{m,k} \tilde{y}_{m,k}}{\sqrt{\sum_{m=1}^{M} \sum_{k=1}^{K_m} (\tilde{x}_{m,k})^2} \sqrt{\sum_{m=1}^{M} \sum_{k=1}^{K_m} (\tilde{y}_{m,k})^2}}$$


### 3. Rare-Variant Generalized Jaccard Similarity

When restricted to rare variants, the Generalized Jaccard formula becomes an incredibly intuitive metric: it calculates the number of *shared rare alleles* divided by the *total unique rare alleles* present across both individuals.

**Mathematical Definition:**
$$J_{rare}(X, Y) = \frac{\sum_{m=1}^{M} \sum_{k=1}^{K_m} \min(\tilde{x}_{m,k}, \tilde{y}_{m,k})}{\sum_{m=1}^{M} \sum_{k=1}^{K_m} \max(\tilde{x}_{m,k}, \tilde{y}_{m,k})}$$


### 4. Mathematical Definition of the Rare-Variant GRM (RV-GRM)


To construct the RV-GRM, we calculate the variance-standardized covariance *only* across the unmasked alleles. The genetic relatedness between individuals $X$ and $Y$ based on rare variants ($A_{XY}^{rare}$) is:

$$A_{XY}^{rare} = \frac{\sum_{m=1}^{M} \sum_{k=1}^{K_m} w_{m,k} (x_{m,k} - 2p_{m,k})(y_{m,k} - 2p_{m,k})}{\sum_{m=1}^{M} \sum_{k=1}^{K_m} w_{m,k} 2p_{m,k}(1 - p_{m,k})}$$
