# -------------------------------------------------------------------------------------------------------------------------------------------------------------
# Script with functions to import and transform raw data
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

"""
    import_data(dir::String, revision::String, year::Integer, N::Integer, S::Integer)

The function loads the WIOD Input-Output Table of specified year as DataFrame into the environment and 
    extracts the raw intermediate, final demand matrices (Z, F) and the inventory adjustment matrix IV.

# Arguments
- `dir::String`: directory of raw data (either folder location or entire directory).
- `revision::String`: specifies the WIOD revision (2013 or 2016) as a string.
- `year::Integer`: specifies the year of the WIOD table which should be imported.
- `N::Integer`: number of origin/destination countries.
- `S::Integer`: number of origin/destination industries.

# Output
- `Z::Matrix`: NS×NS, raw origin country-industry destination country-industry intermediate demand matrix Z.
- `F::Matrix`: NS×N, raw origin country-industry destination country final demand matrix F.
- `IV::Matrix`: NS×N, raw origin country-industry destination country inventory adjustment matrix IV.

# Examples
```julia-repl
julia> Z, F, IV = import_data("C:/Users/u0148308/Desktop/raw/", "2016", 2014, 44, 56)
```
"""

