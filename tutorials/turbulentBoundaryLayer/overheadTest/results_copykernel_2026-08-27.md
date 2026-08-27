# Copy-kernel A/B — A6000 (istmcetus GPU 1), 2026-08-27

Same-day A/B of the same-level halo-copy kernel split (`comm.f90`)
against the pre-change binary `bcc69eb`, both default flags, 400 cold
steps. Generated with `SUFFIX=refbin BIN=<ref>/moby_solve ./run_overhead.sh`
and `./summarise.py`.

| config | reference bcc69eb | candidate | gain |
|---|---|---|---|
| base_jacobi (nb unset) | 1.1836 | 1.1764 | 0.6 % |
| nb16_jacobi | 1.6923 | 1.6118 | 4.8 % |
| refined_yp100_jacobi | 0.7809 | 0.7495 | 4.0 % |
| **block tax (nb16/base)** | **1.4297** | **1.3701** | 14 % of the tax |

The nb-unset row is the control: with one block per rank there are
almost no block-pair halo copies, so there is almost nothing to gain.

```
run                       chron s/step   marginal   ratio    cells   leaves        L2_div     Linf_vel
--------------------------------------------------------------------------------------------------------
base_jacobi                     1.1764     1.1772                         -  1.844810e-05 1.003416e+00
base_jacobi_pre                 1.2422     1.2417                         -  1.844810e-05 1.003416e+00
base_jacobi_prof                1.2412     1.2403                         -  1.844810e-05 1.003416e+00
base_jacobi_refbin              1.1836     1.1845                         -  1.844810e-05 1.003416e+00
base_redblack                   1.4499     1.4529                         -  6.555892e-06 1.003634e+00
base_redblack_prof              1.5083     1.5292                         -  6.555892e-06 1.003634e+00
nb16_jacobi                     1.6118     1.6107   1.368    1.000        -  1.844810e-05 1.003416e+00
nb16_jacobi_pre                 1.7551     1.7533   1.412    1.000        -  1.844810e-05 1.003416e+00
nb16_jacobi_prof                1.7545     1.7522   1.413    1.000        -  1.844810e-05 1.003416e+00
nb16_jacobi_refbin              1.6923     1.6922   1.429    1.000        -  1.844810e-05 1.003416e+00
nb16_redblack                   2.0490     2.0540   1.414    1.000        -  6.555892e-06 1.003634e+00
nb16_redblack_prof              2.0482     2.0535   1.343    1.000        -  6.555892e-06 1.003634e+00
nb8_jacobi_prof                 2.4052     2.4029   1.937    1.000        -  1.844810e-05 1.003416e+00
refined_yp100_jacobi            0.7495     0.7492   0.465    0.455    15360  2.920907e-05 1.002471e+00
refined_yp100_jacobi_pre        0.8670     0.8713   0.497    0.455    15360  2.920907e-05 1.002471e+00
refined_yp100_jacobi_prof        0.8666     0.8709   0.497    0.455    15360  2.920907e-05 1.002471e+00
refined_yp100_jacobi_refbin        0.7809     0.7812   0.462    0.455    15360  2.920907e-05 1.002471e+00

nb-independence nb16_jacobi vs base_jacobi: OK (runtime lines identical)
nb-independence nb16_jacobi_pre vs base_jacobi_pre: OK (runtime lines identical)
nb-independence nb16_jacobi_prof vs base_jacobi_prof: OK (runtime lines identical)
nb-independence nb16_jacobi_refbin vs base_jacobi_refbin: OK (runtime lines identical)
nb-independence nb16_redblack vs base_redblack: OK (runtime lines identical)
nb-independence nb16_redblack_prof vs base_redblack_prof: OK (runtime lines identical)
nb-independence nb8_jacobi_prof vs base_jacobi_prof: OK (runtime lines identical)
```
