# PyPrism

A deconvolution framework for bulk RNA-seq count data based on BayesPrism, implementing the improved algorithm of 
InstaPrism in Python, integrated into the ScVerse, and interoperable with AnnData and HDF5 files.

## Deconvolution

```{python}
deconvolution_result = multi_deconvolution(bulk_data=bulk_data, 
                                           single_cell_reference=reference, 
                                           number_of_iterations=100)
```

