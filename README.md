Simulation studies for Multiset sparse redundancy analysis
========================

Multiset sparse Redundancy Analysis (multi‐sRDA) for multi-source high‐dimensional omics data.

Redundancy Analysis (RDA) is a well‐known method used to describe the directional relationship between related data sets. Recently, we proposed sparse Redundancy Analysis (sRDA) for high‐dimensional genomic data analysis to find explanatory variables that explain the most variance of the response variables. As more and more biomolecular data become available from different biological levels, such as genotypic and phenotypic data from different omics domains, a natural research direction is to apply an integrated analysis approach in order to explore the underlying biological mechanism of certain phenotypes of the given organism. We show that the multiset sparse Redundancy Analysis (multi‐sRDA) framework is a prominent candidate for high‐dimensional omics data analysis since it accounts for the directional information transfer between omics sets, and, through its sparse solutions, the interpretability of the result is improved. In this paper, we also describe a software implementation for multi‐sRDA, based on the Partial Least Squares Path Modeling algorithm. We test our method through simulation and real omics data analysis with data sets of 364,134 methylation markers, 18,424 gene expression markers, and 47 cytokine markers measured on 37 patients with Marfan syndrome.

For more details see [Csala, A., Hof, M. H., & Zwinderman, A. H. (2019). Multiset sparse redundancy analysis for high‐dimensional omics data. Biometrical Journal, 61(2), 406-423.](https://doi.org/10.1002/bimj.201700248).

Please run the code in analyse_multi_sRDA.R to reprudice the results in the manusript. For questions, comments or remarks about the code please contact [a at csala.me].

The code is based on PLS-PM as described in Esposito Vinzi and Russolillo, 2013; Sanchez, 2013.

List of confiurations on which the code was tested:

	R version 3.4.2 (2017-09-28)
	Platform: x86_64-w64-mingw32/x64 (64-bit)
	Running under: Windows 7 x64 (build 7600)

	attached base packages:
	stats     
	graphics  
	grDevices utils     
	datasets  
	methods   
	base     

	other attached packages:
	amap_0.8-14    
	diagram_1.6.4  
	shape_1.4.4    
	turner_0.1.7   
	tester_0.1.7   
	elasticnet_1.1 
	lars_1.2
	plotrix_3.7 


Reference:

	-	Esposito Vinzi, V. and Russolillo, G. (2013). Partial least squares algorithms and methods. Wiley Interdisciplinary
	-	Reviews: Computational Statistics, 5(1):1ñ19.
	Sanchez, G. (2013). Pls path modeling with r. Berkeley: Trowchez Editions.