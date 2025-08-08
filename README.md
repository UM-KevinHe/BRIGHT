BRIGHT
=======================

# Description

`BRIGHT` is a group of methods for using either published genotype-trait summary statistics (GWAS) or individual-level data from different ethnic populations to improve the recalibration, discrimination, and prediction accuracy on the target minority cohort. We implement (group) LASSO, Elastic Net, MCP and SCAD penalties for marker fine mapping.

# Models:
<table>
    <tr>
        <th>Method</th>
        <th>Description</th>
        <th>Example</th>
    </tr>
    <tr>
        <td>BRIGHTs</td>
        <td>
        BRIGHTs utilize summary statistics from both target minority and prior majority data to construct trans-ethnic PRS. Both binary and continuous outcomes are implemented in the package. <a href="#references">[1]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/BRIGHT/">tutorial</a></td>
    </tr>
      <tr>
        <td>BRIGHTi</td>
        <td>
        BRIGHTi utilize individual-level data from target minority population and summery-level data from prior majority population to construct trans-ethnic PRS. Binary, continuous, and count outcomes are implemented in the package. <a href="#references">[2]</a>.
        </td>
        <td><a href="https://um-kevinhe.github.io/BRIGHT/">tutorial</a></td>
    </tr>
</table>

We recommand to start with <a href="https://um-kevinhe.github.io/BRIGHT/" target="_blank">tutorial</a>, where we provide an overview of the package's usage, including preprocessing, model training, hyperparameter fine tunning, and post-estimation evaluation.

# Reference
Li Q, Patrick MT, Zhang H, Khunsriraksakul C, Stuart PE, Gudjonsson JE, Nair R, Elder JT, Liu DJ, Kang J, Tsoi LC, He K. Bregman Divergence-Based Data Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment. arXiv preprint arXiv:2210.06025. 2022 Oct 12. (https://arxiv.org/abs/2210.06025)

Li Q, Patrick MT, Fritsche L, Wang D, Kang J, Tsoi LC, He K. Bregman Divergence-Based Data Integration with Application to Polygenic Risk Score (PRS) Heterogeneity Adjustment. arXiv preprint arXiv:2210.06025. 2022 Oct 12. (https://arxiv.org/abs/2210.06025)