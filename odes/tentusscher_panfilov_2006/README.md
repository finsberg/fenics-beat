# Tentusscher, Panfilov, 2006

To reproduce files in this folder, first get the model files from CellML
```
git clone https://models.physiomeproject.org/workspace/604
```
Install `gotran`
```
python3 -m pip install gotran
```

Convert cellml files to to gotran ode files
```
cellml2gotran ten_tusscher_model_2006_endo.cellml
cellml2gotran ten_tusscher_model_2006_epi.cellml
cellml2gotran ten_tusscher_model_2006_m.cellml
```
This will generate the following files
```
tentusscher_panfilov_2006_endo_cell.ode
tentusscher_panfilov_2006_epi_cell.ode
tentusscher_panfilov_2006_M_cell.ode
```
Now convert them to python code
```¨
gotran2py tentusscher_panfilov_2006_epi_cell.ode --solvers.explicit_euler.generate=1 --solvers.generalized_rush_larsen.generate=1 --code.body.use_enum=1 --namespace=np
gotran2py tentusscher_panfilov_2006_endo_cell.ode --solvers.explicit_euler.generate=1 --solvers.generalized_rush_larsen.generate=1 --code.body.use_enum=1 --namespace=np
gotran2py tentusscher_panfilov_2006_M_cell.ode --solvers.explicit_euler.generate=1 --solvers.generalized_rush_larsen.generate=1 --code.body.use_enum=1 --namespace=np
```
