# Single-cell-RNA-analysis

Series of scripts used for the analysis of single-cell RNA data.
1. SoupX is used to control for ambient RNA, per library.
2. ScDblFinder is used to identify potential doublets.
3. Cell baroced are renamed based on the library they came from. In the same script a merged object is created.
4. Quality control. Cells of low quality are removed either using hard thresholds, or using quantiles per library.
5. Data are normalized, log-transformed, and PCA is performed.
