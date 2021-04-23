File introduction:

1. Parallel_Wave_Equation.cpp: 
the parallel version of Serial_Wave_Equation.cpp

2. out_parallel/:
the files for all cores for each output

3. out_parallel_combined:
the files combined all cores for each output

4. visualization_save.ipynb: 
visualization of results.(You should run the cell again when open the notebook to see the inline animation).
It collectecd the files from out_parallel and combined, then write to out_parallel_combined.
Then it use the files in out_parallel_combined to generate animation in notebook and also saved in current folder.

5. logs/: 
the results running the program on hpc

6. report: 
the report of this coursework. It embeds the visualization code((You should run the cell again when open the notebook to see the inline animation)).

7. img/: 
the images for report

8. dynamic image.gif: 
the animation file


Getting start:

1. Configure the value in Parallel_Wave_Equation.cpp
/*
	Problem config value
*/
//entire gird size
int imax = 100, jmax = 100;

//tmax for iteration
double t_max = 30.0;

//output time index，output time step size
t_out = 0.0, dt_out = 0.04;

//domain size
double y_max = 10.0, x_max = 10.0

 // c value in equation
double c = 1;

2. mpic++ Parallel_Wave_Equation.cpp

3. mpirun -n core_number ./a.out