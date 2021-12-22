# IOTables

[![Build Status](https://github.com/forsthuber92/IOTables.jl/workflows/CI/badge.svg)](https://github.com/forsthuber92/IOTables.jl/actions)
[![Coverage](https://codecov.io/gh/forsthuber92/IOTables.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/forsthuber92/IOTables.jl)


# Scripts

A brief overview about the script and their functions in this repository.

## model

- **transform_data.jl:** transforms raw data used in *baseline_working.jl*
- **price_hat.jl:** inner loop of *baseline_working.jl* to compute price indices
- **wage_hat.jl:** computes wage adjustments by balancing trade for *baseline_working.jl*
- **baseline_working.jl:** summary script to run the baseline model (with or without trade balance)

## counterfactual

- **tariffs_functions.jl:** transforms EU MFN tariffs to be used as counterfactual in exit of EU scenario
- **head_ries_index.jl:** computes symmetric bilateral trade costs according to Head and Ries (2001

## elasticity

- **tariff_data.jl:** transforms raw MFN tariffs for 2018 to be used in *elasticities_functions.jl*
- **elasticities_functions.jl:** computes the statistics in Caliendo and Parro (2015) used in *trade_elasticities.jl*
- **trade_elasticities.jl:** estimates trade elasticity according to Caliendo and Parro (2015)
- **trade_elasticities_EU.jl:** estimates trade elasticity according to Caliendo and Parro (2015) but taking EU as aggregate
