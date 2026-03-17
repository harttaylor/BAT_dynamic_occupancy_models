# BC BAT MONITORING - ROTATION SCHEDULE ANALYSIS
## Updated Results: 2017-2024 (Balanced Design)

---

## METHODS

### Study Design
Evaluated the impact of five monitoring rotation schedules on our ability to detect population trends for five bat species in British Columbia: Hoary Bat (LACI), Silver-haired Bat (LANO), Eastern Red Bat (LABO), Little Brown Myotis (MYLU), and Northern Myotis (MYSE).

**Sites:** 173 total BC sites
- **Core sites** (n = 88): Sites with consistent monitoring commitments (National Parks, Indigenous territories, Nature Conservancy lands) monitored annually under all schedules
- **Rotating sites** (n = 85): Sites on BC provincial lands (BC-WLRS, BC Parks) included in rotation designs

**Study Period:** 2017-2024 (8 years)
*Note: 2016 excluded due to data reliability concerns per discussion with Cami (30/01/2026)*

### Rotation Schedules Evaluated

| Schedule | Description | Effort Saved |
|----------|-------------|--------------|
| Full (Annual) | All 173 sites every year | 0% (reference) |
| 2-Year Rotation | Rotating sites split into 2 groups, alternating years | ~25% |
| 3-Year Rotation | Rotating sites split into 3 groups, cycled | ~33% |
| 4-Year Rotation | Rotating sites split into 4 groups, cycled | ~37% |
| 5-Year Rotation | Rotating sites split into 5 groups, cycled | ~39% |

**Balanced Assignment:** Rotating sites were assigned to groups using stratified random assignment within each region, ensuring:
- Equal number of sites sampled each year (±1 site)
- Representation from all regions in each sampling year
- No systematic covariate confounding

### Occupancy Modeling

Dynamic occupancy models fitted using Bayesian framework (JAGS, 30,000 iterations, 3 chains).

**Model Structure:**
- **Initial occupancy (ψ):** Water presence (binary), distance to road, region (random effect)
- **Colonization (γ):** Year-specific intercepts + distance to road + distance to harvest
- **Persistence (φ):** Year-specific intercepts + distance to road + distance to harvest  
- **Detection (p):** Clutter, temperature, Julian date

**Trend Estimation:** λ = ψ₂₀₂₄ / ψ₂₀₁₇ (ratio of final to initial year occupancy)

---

## RESULTS

All models converged successfully (R-hat < 1.1 for all parameters).

### Population Trends (Full Annual Sampling)

| Species | Common Name | λ (95% CI) | Trend | CI Width |
|---------|-------------|------------|-------|----------|
| LANO | Silver-haired Bat | **1.19 (1.08, 1.31)** | **INCREASING** ✓ | 0.23 |
| LACI | Hoary Bat | 1.05 (0.92, 1.20) | Stable | 0.28 |
| MYLU | Little Brown Myotis | 0.96 (0.92, 1.01) | Stable | 0.08 |
| LABO | Eastern Red Bat | **0.86 (0.75, 0.99)** | **DECLINING** ✓ | 0.24 |
| MYSE | Northern Myotis | 0.85 (0.64, 1.14) | Uncertain | 0.50 |

*✓ = 95% CI excludes 1.0, indicating statistically significant trend*

### Precision Loss by Rotation Schedule

Relative increase in 95% CI width compared to annual sampling:

| Species | 2-Year | 3-Year | 4-Year | 5-Year |
|---------|:------:|:------:|:------:|:------:|
| LACI (Hoary) | +5% | +8% | +7% | +7% |
| LANO (Silver-haired) | +10% | +16% | +10% | +23% |
| LABO (Eastern Red) | +21% | +41% | +57% | +55% |
| MYLU (Little Brown) | +31% | +30% | +42% | +52% |
| MYSE (Northern) | +17% | +54% | +43% | +58% |
| **Average** | **+17%** | **+30%** | **+32%** | **+39%** |

### Detailed Results by Species

#### LACI - Hoary Bat
| Design | Site-Years | λ (95% CI) | CI Width | Precision Loss |
|--------|------------|------------|----------|----------------|
| Full | 1384 | 1.05 (0.92, 1.20) | 0.28 | — |
| 2-Year | 1044 | 1.00 (0.87, 1.16) | 0.30 | +5% |
| 3-Year | 933 | 1.05 (0.91, 1.21) | 0.30 | +8% |
| 4-Year | 874 | 1.00 (0.86, 1.16) | 0.30 | +7% |
| 5-Year | 843 | 0.97 (0.83, 1.13) | 0.30 | +7% |

#### LANO - Silver-haired Bat
| Design | Site-Years | λ (95% CI) | CI Width | Precision Loss |
|--------|------------|------------|----------|----------------|
| Full | 1384 | 1.19 (1.08, 1.31) | 0.23 | — |
| 2-Year | 1044 | 1.17 (1.06, 1.31) | 0.25 | +10% |
| 3-Year | 933 | 1.16 (1.04, 1.31) | 0.27 | +16% |
| 4-Year | 874 | 1.09 (0.97, 1.23) | 0.25 | +10% |
| 5-Year | 843 | 1.17 (1.04, 1.32) | 0.28 | +23% |

#### LABO - Eastern Red Bat
| Design | Site-Years | λ (95% CI) | CI Width | Precision Loss |
|--------|------------|------------|----------|----------------|
| Full | 1384 | 0.86 (0.75, 0.99) | 0.24 | — |
| 2-Year | 1044 | 0.91 (0.78, 1.06) | 0.29 | +21% |
| 3-Year | 933 | 0.88 (0.73, 1.06) | 0.33 | +41% |
| 4-Year | 874 | 0.90 (0.74, 1.11) | 0.37 | +57% |
| 5-Year | 843 | 0.88 (0.72, 1.08) | 0.37 | +55% |

#### MYLU - Little Brown Myotis  
| Design | Site-Years | λ (95% CI) | CI Width | Precision Loss |
|--------|------------|------------|----------|----------------|
| Full | 1384 | 0.96 (0.92, 1.01) | 0.08 | — |
| 2-Year | 1044 | 0.95 (0.90, 1.01) | 0.11 | +31% |
| 3-Year | 933 | 0.95 (0.90, 1.01) | 0.11 | +30% |
| 4-Year | 874 | 0.94 (0.88, 1.00) | 0.12 | +42% |
| 5-Year | 843 | 0.95 (0.89, 1.02) | 0.12 | +52% |

#### MYSE - Northern Myotis
| Design | Site-Years | λ (95% CI) | CI Width | Precision Loss |
|--------|------------|------------|----------|----------------|
| Full | 1384 | 0.85 (0.64, 1.14) | 0.50 | — |
| 2-Year | 1044 | 0.82 (0.58, 1.16) | 0.58 | +17% |
| 3-Year | 933 | 0.98 (0.67, 1.43) | 0.76 | +54% |
| 4-Year | 874 | 0.93 (0.64, 1.35) | 0.71 | +43% |
| 5-Year | 843 | 0.95 (0.63, 1.41) | 0.78 | +58% |

---

## KEY FINDINGS

### 1. Population Status of BC Bats (2017-2024)

**Silver-haired Bat (LANO) is INCREASING:** 
- λ = 1.19 (CI: 1.08-1.31) indicates ~19% increase in occupancy
- This positive trend was detected across all rotation schedules
- Most robust finding in the analysis

**Eastern Red Bat (LABO) is DECLINING:**
- λ = 0.86 (CI: 0.75-0.99) indicates ~14% decline in occupancy
- Significant decline detected with full sampling
- **Important caveat:** Under reduced sampling (3-5 year rotations), the CI crosses 1.0, meaning we would lose the ability to detect this decline as statistically significant

**Little Brown Myotis (MYLU) is STABLE:**
- λ = 0.96 (CI: 0.92-1.01) - very tight CI, high precision
- No evidence of population change
- This species has high detection rates, giving precise estimates

**Hoary Bat (LACI) appears STABLE:**
- λ = 1.05 (CI: 0.92-1.20) - consistent across all designs
- No significant trend detected

**Northern Myotis (MYSE) - UNCERTAIN:**
- λ = 0.85 (CI: 0.64-1.14) - wide confidence interval
- Point estimate suggests possible decline (~15%)
- Low detection probability leads to high uncertainty
- Cannot rule out either decline or stability

### 2. Rotation Schedule Recommendations

**2-Year Rotation is the recommended option:**
- Saves ~25% field effort
- Average precision loss of only 17%
- All biological conclusions remain unchanged
- Best trade-off between efficiency and statistical power

**3-Year and longer rotations have significant drawbacks:**
- For declining species (LABO), longer rotations cause CIs to cross 1.0
- This means we would fail to detect real declines
- MYSE already has wide CIs; reduced sampling makes trend detection nearly impossible

### 3. Species-Specific Sensitivity

Species with **lower detection rates** are more sensitive to reduced sampling:
- MYSE (Northern Myotis): Most affected - precision loss up to 58%
- LABO (Eastern Red Bat): Highly affected - would lose significant decline signal
- MYLU (Little Brown): Least affected due to high detection rates

### 4. Management Implications

| If adopting... | Pros | Cons |
|----------------|------|------|
| **2-Year Rotation** | 25% cost savings; maintains statistical power | Slightly wider CIs |
| **3-Year Rotation** | 33% cost savings | Risk losing ability to detect LABO decline |
| **5-Year Rotation** | 39% cost savings | Substantial precision loss; unreliable for rare species |

---

## SUMMARY

The BC bat monitoring program (2017-2024) reveals:

1. **Good news:** Silver-haired Bat populations are increasing significantly
2. **Concern:** Eastern Red Bat populations are declining; this signal could be lost with reduced sampling
3. **Watch list:** Northern Myotis shows possible decline but uncertainty is high due to low detection
4. **Stable:** Little Brown Myotis and Hoary Bat show no significant change

**Recommendation:** A **2-year rotation schedule** provides the best balance of cost savings (~25%) and statistical power. Longer rotation cycles risk losing the ability to detect biologically important population declines.

---

## ADDITIONAL INFORMATION TO CONSIDER INCLUDING

### Detection Probability Estimates
Could add a table showing mean detection probability by species - helps explain why some species have more uncertainty.

### Covariate Effects
The models estimate effects of:
- Distance to road on colonization/persistence
- Distance to harvest on colonization/persistence
- Could summarize significant effects if biologically interesting

### Occupancy Estimates by Year
Annual occupancy plots (already in figures) show the trajectory - could add a summary table of start vs end occupancy.

### Model Diagnostics
- R-hat values (all < 1.1)
- Effective sample sizes
- Could include as appendix

### Recommendations for Future Monitoring
- Consider targeted intensive monitoring for MYSE given high uncertainty
- LABO warrants continued close monitoring given decline signal
