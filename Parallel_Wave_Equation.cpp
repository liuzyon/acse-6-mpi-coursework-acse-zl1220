#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

int id, p;
clock_t mpi_start;
clock_t mpi_end;

/*
	This method is used to calculate the rows and columns of processes. 
	For example, if p is 4, it will make rows=2, columns=2 . It arranges 
	the process as close to the square as possible.
*/
void find_dimensions(int p, int &rows, int &columns)
{
	int min_gap = p;
	int top = sqrt(p) + 1;
	for (int i = 1; i <= top; i++)
	{
		if (p % i == 0)
		{
			int gap = abs(p / i - i);

			if (gap < min_gap)
			{
				min_gap = gap;
				rows = i;
				columns = p / i;
			}
		}
	}

	if (id == 0)
	{
		cout << "The total number of processor is: " << p << "." << endl;
		cout << "Divide " << p << " into " << rows << " by " << columns << " grid" << endl;
		mpi_start = clock();
	}
}

// row and column numbers for all processes 
int rows, columns;

// id of row and column for current process
int id_row, id_column;

/*
	This method is used to get the row id and column id by the process id.
*/
void id_to_index(int id, int &id_row, int &id_column)
{
	id_column = id % columns;
	id_row = id / columns;
}

/*
	This method is used to get the process id by the row id and column id.
*/
int id_from_index(int id_row, int id_column)
{
	if (id_row >= rows || id_row<0)
		return -1;
	if (id_column >= columns || id_column<0)
		return -1;

	return id_row*columns + id_column;
}

/*
	This method is used to allocate the subgrid to each process. 
	The rows and columns of each subgird is average value of the 
	whole domain with extra two rows and two columns(ghost layers 
	for communication). If there are remaining rows and columns 
	when the whole domain is equally allocated to each process, 
	the remaining row and columns are allocated to rightmost and 
	bottommost side processes.
*/
void allocate_grid(vector< vector<double> > &grid, vector< vector<double> > &new_grid, vector< vector<double> > &old_grid, int imax, int jmax, int &grid_rows, int &grid_columns)
{
	grid_rows = imax / rows + 2;
	grid_columns = jmax / columns + 2;
	int rows_remain = imax % rows;
	int columns_remain = jmax % columns;

	if (id_row == rows-1 && id_column == columns-1)
	{
		// The process in the last row and last column allocates the remaining rows and columns
		old_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns+columns_remain));
		grid.resize(grid_rows+rows_remain, vector<double>(grid_columns+columns_remain));
		new_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns+columns_remain));

	} else if (id_row == rows-1)
	{	// The process on the last row allocates the remaining rows
		old_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns));
		grid.resize(grid_rows+rows_remain, vector<double>(grid_columns));
		new_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns));

	} else if (id_column == columns-1)
	{
		// The process on the last column allocates the remaining columns
		old_grid.resize(grid_rows, vector<double>(grid_columns+columns_remain));
		grid.resize(grid_rows, vector<double>(grid_columns+columns_remain));
		new_grid.resize(grid_rows, vector<double>(grid_columns+columns_remain));
	} else
	{	// Other processes are assigned on an average row and column
		old_grid.resize(grid_rows, vector<double>(grid_columns));
		grid.resize(grid_rows, vector<double>(grid_columns));
		new_grid.resize(grid_rows, vector<double>(grid_columns));
	}
	
}

/*
	This method is used to sets half sinusoidal intitiak disturbance 
	for subdomain for its corresponding position in the entire domain
*/
void init(vector< vector<double> > &grid, vector< vector<double> > &old_grid, double dx, double dy, int id_row, int id_column, int grid_rows, int grid_columns)
{
	//sets half sinusoidal intitial disturbance
	//average gird row and column (exclude ghost layer)
	int avg_rows = grid_rows-2;
	int avg_columns = grid_columns-2;

	//current grid row and column(include ghost layer)
	int p_grid_rows = grid.size();
	int p_grid_columns = grid[0].size();

	double r_splash = 1.0;
	double x_splash = 3.0;
	double y_splash = 3.0;

	for (int i = 1; i < p_grid_rows - 1; i++)
		for (int j = 1; j < p_grid_columns - 1; j++)
		{
			double x = dx * (i-1) + dx * avg_rows * id_row;
			double y = dy * (j-1) + dy * avg_columns * id_column;

			double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));

			if (dist < r_splash)
			{
				double h = 5.0*(cos(dist / r_splash * M_PI) + 1.0);

				grid[i][j] = h;
				old_grid[i][j] = h;
			}
		}

}

/*
	Write the gird value to output file for this process, not write the ghost layer.
*/
void grid_to_file(int id, int out, vector< vector<double> > &grid)
{
	int imax = grid.size();
	int jmax = grid[0].size();
	
	stringstream fname;
	fstream f1;
	fname << "./out_parallel/output" << "_" << out << "_p" << id << ".dat";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < imax-1; i++)
	{
		for (int j = 1; j < jmax-1; j++)
			f1 << grid[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}

// MPI_Datatype for sending
MPI_Datatype Datatype_left, Datatype_right, Datatype_top, Datatype_bottom;

/*
	Create the MPI data type for sending, which is Datatype_left, Datatype_right, 
	Datatype_top, Datatype_bottom. The four data type correspond to the four edge 
	of data_array.
*/
void createdatatypes(double** data, int m, int n)
{
	vector<int> block_lengths;
	vector<MPI_Datatype> typelist;
	vector<MPI_Aint> addresses;
	MPI_Aint add_start;

	//left
	for (int i = 0; i < m; i++)
	{
		block_lengths.push_back(1);
		typelist.push_back(MPI_DOUBLE);
		MPI_Aint temp_address;
		MPI_Get_address(&data[i][0], &temp_address);
		addresses.push_back(temp_address);
	}
	MPI_Get_address(data, &add_start);
	for (int i = 0; i < m; i++) addresses[i] = addresses[i] - add_start;
	MPI_Type_create_struct(m, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_left);
	MPI_Type_commit(&Datatype_left);

	//right
	block_lengths.resize(0);
	typelist.resize(0);
	addresses.resize(0);
	for (int i = 0; i < m; i++)
	{
		block_lengths.push_back(1);
		typelist.push_back(MPI_DOUBLE);
		MPI_Aint temp_address;
		MPI_Get_address(&data[i][n - 1], &temp_address);
		addresses.push_back(temp_address);
	}
	for (int i = 0; i < m; i++) addresses[i] = addresses[i] - add_start;
	MPI_Type_create_struct(m, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_right);
	MPI_Type_commit(&Datatype_right);

	//top - only need one value
	int block_length = n;
	MPI_Datatype typeval = MPI_DOUBLE;
	MPI_Aint address;
	MPI_Get_address(data[0], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_top);
	MPI_Type_commit(&Datatype_top);

	//bottom - only need one value
	MPI_Get_address(data[m - 1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_bottom);
	MPI_Type_commit(&Datatype_bottom);
}

// MPI_Datatype for receiving
MPI_Datatype Datatype_left_recv, Datatype_right_recv, Datatype_top_recv, Datatype_bottom_recv;

/*
	Create the MPI data type for receiving, which is Datatype_left_recv, 
	Datatype_right_recv, Datatype_top_recv, Datatype_bottom_recv. The four 
	data type correspond to the four edge of recv_array
*/
void createdatatypes_recv(double** data, int m, int n)
{
	vector<int> block_lengths;
	vector<MPI_Datatype> typelist;
	vector<MPI_Aint> addresses;
	MPI_Aint add_start;

	// left
	for (int i = 1; i < m-1; i++)
	{
		block_lengths.push_back(1);
		typelist.push_back(MPI_DOUBLE);
		MPI_Aint temp_address;
		MPI_Get_address(&data[i][0], &temp_address);
		addresses.push_back(temp_address);
	}
	MPI_Get_address(data, &add_start);
	for (int i = 0; i < m-2; i++) addresses[i] = addresses[i] - add_start;
	MPI_Type_create_struct(m-2, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_left_recv);
	MPI_Type_commit(&Datatype_left_recv);

	// right
	block_lengths.resize(0);
	typelist.resize(0);
	addresses.resize(0);
	for (int i = 1; i < m-1; i++)
	{
		block_lengths.push_back(1);
		typelist.push_back(MPI_DOUBLE);
		MPI_Aint temp_address;
		MPI_Get_address(&data[i][n - 1], &temp_address);
		addresses.push_back(temp_address);
	}
	for (int i = 0; i < m-2; i++) addresses[i] = addresses[i] - add_start;
	MPI_Type_create_struct(m-2, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_right_recv);
	MPI_Type_commit(&Datatype_right_recv);

	//top - only need one value
	int block_length = n-2;
	MPI_Datatype typeval = MPI_DOUBLE;
	MPI_Aint address;
	MPI_Get_address(&data[0][1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_top_recv);
	MPI_Type_commit(&Datatype_top_recv);

	//bottom - only need one value
	MPI_Get_address(&data[m - 1][1], &address);
	address = address - add_start;
	MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_bottom_recv);
	MPI_Type_commit(&Datatype_bottom_recv);


}

/*
	Copy the data of edges in recv_array 
	to gird(only copy the edge it received, 
	not all four edges).
*/
void copy_boundry(vector< vector<double> > &grid, int recv_array_rows, int recv_array_columns, double** recv_array, bool recv_left, bool recv_top, bool recv_right, bool recv_bottom)
{
	// gird size
	int imax = grid.size();
	int jmax = grid[0].size();

	if (recv_left)
	{	
		// left ghost layer
		for (int i = 1; i < imax-1; i++)
		{
			grid[i][0] = recv_array[i][0];
		}
	}
	
	if (recv_right)
	{
		// right ghost layer
		for (int i = 1; i < imax-1; i++)
		{
			grid[i][jmax-1] = recv_array[i][recv_array_columns-1]; // 这里需要改
		}
	}
	
	if (recv_top)
	{
		// top ghost layer
		for (int i = 1; i < jmax-1; i++)
		{
			grid[0][i] = recv_array[0][i];
		}
	}
	

	if (recv_bottom)
	{
		// bottom ghost layer
		for (int i = 1; i < jmax-1; i++)
		{
			grid[imax-1][i] = recv_array[recv_array_rows-1][i];
		}
	}

}

/*
	Problem config value
*/
//entire gird size
int imax = 100, jmax = 100;

//gird for each process
vector< vector<double> > grid, new_grid, old_grid;

//tmax for iteration
double t_max = 30.0;

//time index，time step size，output time index，output time step size
double t, dt, t_out = 0.0, dt_out = 0.04;

//domain size
double y_max = 10.0, x_max = 10.0, dx, dy;

 // c value in equation
double c = 1;

/* 
	Perform the calculation for new_gird using grid and old_grid, 
	and set the boundry value(Neumann boundaries) if process is in boundry.
*/
void do_iteration()
{
	//gird size (include ghost layer)
	int imax = grid.size();
	int jmax = grid[0].size();

	// if process is not in boundry
	if (id_row != 0 && id_row != rows-1 && id_column != 0 && id_column != columns-1)
	{
		for (int i = 1; i < imax - 1; i++)
			for (int j = 1; j < jmax - 1; j++)
				new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];
	}

	// if process is in top boundry
	if (id_row == 0)
	{
		if (id_column == 0) // if process is in left boundry
		{
			for (int i = 2; i < imax - 1; i++)	// not iterate top boundry
				for (int j = 2; j < jmax - 1; j++)	//not iterate left boundry
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else if (id_column == columns-1) // if process is in right boundry
		{
			for (int i = 2; i < imax - 1; i++)	// not iterate top boundry
				for (int j = 1; j < jmax - 2; j++)	// not iterate right boundry
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else	// if process is both in top boundry and left boundry
		{
			for (int i = 2; i < imax - 1; i++)	// not iterate top boundry
				for (int j = 1; j < jmax - 1; j++)	// iterate left boundry and right boundry
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];
		}

	}

	// if process is not in bottom boundry
	if (id_row == rows-1)
	{
		if (id_column == 0)	// if process is in left boundry
		{
			for (int i = 1; i < imax - 2; i++)	// not iterate bottom boundry
				for (int j = 2; j < jmax - 1; j++)	//not iterate left boundry
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else if (id_column == columns-1) // if process is in right boundry
		{
			for (int i = 1; i < imax - 2; i++)	// not iterate bottom boundry
				for (int j = 1; j < jmax - 2; j++)	// not iterate right boundry
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else	// if process is both in top boundry and left boundry
		{
			for (int i = 1; i < imax - 2; i++)	// not iterate bottom boundry
				for (int j = 1; j < jmax - 1; j++) // iterate left boundry and right boundry
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];
		}

	}
	
	// set Neumann boundaries according process position

	// if it is in top boundry
	if(id_row == 0)
	{
		// set the top Neumann boundaries
		for (int j = 1; j < jmax-1; j++)
		{
			new_grid[1][j] = new_grid[2][j];
		}

	}


	// if it is in bottom boundry
	if(id_row == rows-1)
	{
		// set the bottom Neumann boundaries
		for (int j = 1; j < jmax-1; j++)
		{
			new_grid[imax-2][j] = new_grid[imax-3][j];
		}
	}

	// if it is in left boundry
	if(id_column == 0)
	{
		// set the left Neumann boundaries
		for (int i = 1; i < imax-1; i++)
		{
			new_grid[i][1] = new_grid[i][2];
		}

	}

	// if it is in right boundry
	if(id_column == columns-1)
	{
		// set the right Neumann boundaries
		for (int i = 1; i < imax-1; i++)
		{
			new_grid[i][jmax - 2] = new_grid[i][jmax - 3];
		}
	}

	t += dt;

	old_grid.swap(new_grid);
	old_grid.swap(grid);
}


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);	//process id
	MPI_Comm_size(MPI_COMM_WORLD, &p);	//process number

	// calculate the rows and columns of processes according to p
	find_dimensions(p, rows, columns);

	// get the row id and column id by the process id
	id_to_index(id, id_row, id_column);

	// average grid row and column value(include ghost layer)
	int grid_rows, grid_columns;

	// allocate the subgrid to each process
	allocate_grid(grid, new_grid, old_grid, imax, jmax, grid_rows, grid_columns);

	dx = x_max / ((double)imax - 1);	//delta x
	dy = y_max / ((double)jmax - 1);	//delta y

	t = 0.0;

	dt = 0.1 * min(dx, dy) / c; // delta t

	int out_cnt = 0, it = 0;	//output time index，iterate time index in iteration

	// sets half sinusoidal intitiak disturbance for subdomain for its corresponding position in the entire domain.
	init(grid, old_grid, dx, dy, id_row, id_column, grid_rows, grid_columns);

	// wrtie the initial output
	grid_to_file(id, out_cnt, grid);
	out_cnt++;
	t_out += dt_out;

	//Buffer used to sending and recving
	double** data_array;
	double** recv_array;


	// Build the send and receive array memory

    // Use the current process grid to create the sending array(include ghost layer)
	const int m = grid.size()-2;
	const int n = grid[0].size()-2;
	data_array = new double*[m];
	for (int i = 0; i < m; i++)
		data_array[i] = new double[n];

	// Use the max subgrid size to build receiving array to prevent truncation 
	// error when receiving(exclude ghost layer).
	int rows_remain = imax % rows;
	int columns_remain = jmax % columns;
	const int l = grid_rows+rows_remain;
	const int k = grid_columns+columns_remain;
	recv_array = new double*[l];
	for (int i = 0; i < l; i++)
		recv_array[i] = new double[k];


	// Create the MPI data type for sending
	createdatatypes(data_array, m, n);
	// Create the MPI data type for receiving
	createdatatypes_recv(recv_array, l, k);

	while (t < t_max)
	{
		
		// Ready for send data
        int gr = grid.size() - 2; // current grid valid rows(exclude ghost layer)
        int gc = grid[0].size() - 2;    //current grid valid columns(exclude ghost layer)
		for (int i = 0; i < gr; i++)
			for (int j = 0; j < gc; j++)
				data_array[i][j] =(double) grid[i+1][j+1]; //copy the grid data to sending buffer

		// Initialize the recv_array
		for (int i = 0; i < l; i++)
			for (int j = 0; j < k; j++)
				recv_array[i][j] = -1;

		//sending count
		int cnt = 0;

		// Requsets memory
		MPI_Request* request = new MPI_Request[4 * 2];

		bool recv_top = false;
		bool recv_left = false;
		bool recv_right = false;
		bool recv_bottom = false;
		
		// Iterate the neighbor processes
		for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
			{
				// The process sended row index and column index
				int com_i = id_row + i;
				int com_j = id_column + j;

				// The process sended id
				int com_id = id_from_index(com_i, com_j);

				// If process sended is not current process id, and it is an existed id
				if (com_id != id && com_id >= 0 && com_id < p)
				{
					// If process sended is in the top of current process
					if (com_i == id_row-1 && com_j == id_column)
					{
						// tag_num set to iterate index
						MPI_Isend(data_array, 1, Datatype_top, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_top_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						// recv_top set to true if current process has top process to receive data
						recv_top = true;
						cnt++;
					} 
					else if (com_i == id_row && com_j == id_column-1) // If it is in the left of current process
					{
						MPI_Isend(data_array, 1, Datatype_left, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_left_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						// recv_left set to true if current process has left process to receive data
						recv_left = true;
						cnt++;
					} 
					else if (com_i == id_row && com_j == id_column+1)// If it is in the right of right process
					{
						MPI_Isend(data_array, 1, Datatype_right, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_right_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						// recv_right set to true if current process has right process to receive data
						recv_right = true;
						cnt++;
					}
					else if (com_i == id_row+1 && com_j == id_column) // If it is in the bottom of current process
					{
						MPI_Isend(data_array, 1, Datatype_bottom, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_bottom_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						// recv_bottom set to true if current process has bottom process to receive data
						recv_bottom = true;
						cnt++;
					}

				}
			}

		// Wait for all sending and receiving for this current process finished
		MPI_Waitall(cnt * 2, request, MPI_STATUS_IGNORE);

		// Copy the data of edges in recv_array to gird(only copy the edge it received, not all four edges).
		copy_boundry(grid, l, k, recv_array, recv_left, recv_top, recv_right, recv_bottom);

		// Perform the calculation for new_gird using grid and old_grid, and set the boundry value(Neumann boundaries) if process is in boundry.
		do_iteration();

		// Write the output if it should
		if (t_out <= t)
		{
			grid_to_file(id, out_cnt, grid);
			out_cnt++;
			t_out += dt_out;
		}

		it++;

		delete[] request;
	}

	// Calculate the program run time
	if (id == 0)
	{
		mpi_end = clock();
		double time = double(mpi_end-mpi_start)/CLOCKS_PER_SEC;
		cout << "Total time consumption: " << time << endl;
	}
	

	MPI_Finalize();

	delete[] data_array;
	delete[] recv_array;
}
