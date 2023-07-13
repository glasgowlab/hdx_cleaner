# HDX cleaner and parser

**forked from malcolm's script**

## Main purpose

* clean/parse the raw output from HDExaminer to make uptake plots, pymol plots and etc.

## Useage

```bash
pymol -r main.py --t peptide_pool_20230502.csv --r rangeslist_e.csv --pm 1pfk_Xray_renum.pdb --cbarmax 0.05 
```

where:
> -r: the main script,
> --t: peptide pool from HDExaminer,
> --r: rangelist of the good peptide,
> -pm: pdb file used to make pymol plot,
> --e: exclude the rangelist as bad peptide
