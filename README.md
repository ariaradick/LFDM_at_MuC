# LFDM_at_MuC
Code to reproduce plots for paper [[2312.03826]](https://arxiv.org/abs/2312.03826).

## Steps:

1. Install [Julia](https://julialang.org/), as all files (including `.ipynb` files) use Julia v1.9.3.

2. Open a terminal and navigate to `one_flav_env/` and run `julia --project=.`. This launches Julia with the virtual environment `one_flav_env`. Once Julia opens, press `]` to go into package mode and run `instantiate`. This will install all the packages necessary to run this code.

3. You now need to run all of the madgraph scripts. `madgraph/madgraph_scripts.jl` will automatically create and run such scripts and move the outputs to `madgraph/data/`, but you'll need to make sure to input the location of your Madgraph executable `MADGRAPH_EXE` and the location where Madgraph should save the outputs initially `MADGRAPH_OUTPUT` (this is because Madgraph can't handle outputting to folders with spaces). Then run `/summarize_runs.jl`. Alternatively, you can download all the madgraph data [from dropbox](https://www.dropbox.com/scl/fo/at0eaguxdbs64m5s8pt7r/h?rlkey=7qdl6je8z5ci0c0zquroi4nuh&dl=0) and make sure to put the `data/` folder into `madgraph/`.

4. Now you can run whatever `.ipynb` files you want.