# Generating equivalence classes of the pallet loading problem

This code generates the minimum size instance of every equivalence class of the pallet loading problem with upto _N_ boxes. For further information on the algorithm, see the two papers below. The spreadsheet summary\_data.xls contains supplemental data for the first paper which this code accompanies.


"An algorithm for identifying equivalence classes of the pallet loading problem" (available at https://ssrn.com/abstract=3964980).
 
"Bounds on equivalence classes of the pallet loading problem" (available at https://ssrn.com/abstract=3964352).


## Usage

To compile, 

```
f95 main.f95
```

This implementation uses static arrays and 64bit integers. As uploaded, it can be run for _N_ upto 600. For larger _N_, increase the parameter _N_MAX_ in main.f95 and recompile.

Nov 16, 2021

