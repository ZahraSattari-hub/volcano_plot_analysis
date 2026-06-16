Differential Protein Abundance Analysis and Volcano Plot <br/>

This analysis evaluates differences in protein abundance of secretome samples from hormone treated and untreated mammary epithelial cell cultures (originally published at https://doi.org/10.1016/j.fufo.2024.100395). <br/>

The package used for the differential analysis was Limma (Linear Models for Microarray and RNA-Seq Data). <br/>

Workflow: <br/>
1.	Protein abundance data were imported from the proteomics dataset and sample metadata were used to define treatment groups.
2.	A design matrix describing the experimental conditions was generated.
3.	Protein-wise linear models were fitted using the limma framework.
4.	Empirical Bayes moderation was applied to improve variance estimation across proteins.
5.	Log2 fold changes, p-values, and multiple-testing adjusted p-values were calculated for all detected proteins.
6.	Proteins were classified as upregulated, downregulated, or not significantly changed using predefined fold-change and significance thresholds.
7.	Results were visualized using a volcano plot displaying treatment-associated changes in protein abundance.
8.	Selected proteins of biological interest were annotated and labelled directly on the figure.
   
The workflow demonstrates a standard differential abundance analysis approach commonly used in proteomics studies and provides a reproducible framework for identifying proteins associated with experimental treatments.

[Plot.tiff](https://github.com/user-attachments/files/28990844/Plot.tiff)

