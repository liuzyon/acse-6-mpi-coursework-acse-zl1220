#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

int id, p, tag_num = 1;

void find_dimensions(int p, int &rows, int &columns)		//A bit brute force - this can definitely be made more efficient!
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
		cout << "Divide " << p << " into " << rows << " by " << columns << " grid" << endl;
}

int rows, columns;
int id_row, id_column;

//通过进程的id得到其grid行列index
void id_to_index(int id, int &id_row, int &id_column)
{
	id_column = id % columns;
	id_row = id / columns;
}

int id_from_index(int id_row, int id_column)
{
	if (id_row >= rows || id_row<0)
		return -1;
	if (id_column >= columns || id_column<0)
		return -1;

	return id_row*columns + id_column;	//通过行列index返回id
}

void allocate_grid(vector< vector<double> > &grid, vector< vector<double> > &new_grid, vector< vector<double> > &old_grid, int imax, int jmax, int &grid_rows, int &grid_columns)
{
	//创建棋盘，按grid行列/进程行列划分，不能整除的余数被第一行第一列进程分配，每个进程同时多分配两行两列作为ghost layer

	grid_rows = imax / rows + 2;
	grid_columns = jmax / columns + 2;
	int rows_remain = imax % rows;
	int columns_remain = jmax % columns;

	if (id_row == rows-1 && id_column == columns-1)
	{
		//最后一行最后一列的进程要分配剩余的行和列
		old_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns+columns_remain));
		grid.resize(grid_rows+rows_remain, vector<double>(grid_columns+columns_remain));
		new_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns+columns_remain));

	} else if (id_row == rows-1)
	{	//最后一行的进程要分配上剩余的行
		old_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns));
		grid.resize(grid_rows+rows_remain, vector<double>(grid_columns));
		new_grid.resize(grid_rows+rows_remain, vector<double>(grid_columns));

	} else if (id_column == columns-1)
	{
		//最后一列的进程要分配上剩余的列
		old_grid.resize(grid_rows, vector<double>(grid_columns+columns_remain));
		grid.resize(grid_rows, vector<double>(grid_columns+columns_remain));
		new_grid.resize(grid_rows, vector<double>(grid_columns+columns_remain));
	} else
	{	//其他进程按平均值分配
		old_grid.resize(grid_rows, vector<double>(grid_columns));
		grid.resize(grid_rows, vector<double>(grid_columns));
		new_grid.resize(grid_rows, vector<double>(grid_columns));
	}
	
}

void init(vector< vector<double> > &grid, vector< vector<double> > &old_grid, double dx, double dy, int id_row, int id_column, int grid_rows, int grid_columns)
{
	//sets half sinusoidal intitial disturbance （设置半正弦波初始扰动） - this is brute force - it can be done more elegantly
	//gird 平均行列数(去除ghost layer)
	int avg_rows = grid_rows-2;
	int avg_columns = grid_columns-2;

	// cout << "dx: " << dx << "dy: " << dy << endl;

	//当前进程grid的行数和列数(包括了ghost layer)
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

MPI_Datatype Datatype_left, Datatype_right, Datatype_top, Datatype_bottom;
// MPI_Datatype Datatype_left;

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

MPI_Datatype Datatype_left_recv, Datatype_right_recv, Datatype_top_recv, Datatype_bottom_recv;

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

void copy_boundry(vector< vector<double> > &grid, int recv_array_rows, int recv_array_columns, double** recv_array, bool recv_left, bool recv_top, bool recv_right, bool recv_bottom)
{
	int imax = grid.size(); // gird 规模
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

int imax = 12, jmax = 12;	//棋盘大小
//进程负责的grid
vector< vector<double> > grid, new_grid, old_grid;
double t_max = 30.0;	//迭代t的最大值
double t, dt, t_out = 0.0, dt_out = 0.04;	//时间index，时间间隔，输出时间index，输出时间间隔 0.04
double y_max = 10.0, x_max = 10.0, dx, dy;
double c = 1;



//进行一次迭代
void do_iteration()
{
	//包含ghost layer
	int imax = grid.size();
	int jmax = grid[0].size();

	// cout << "p" << id << id_row << id_column << rows << columns << endl;

	// 如果不为边界进程
	if (id_row != 0 && id_row != rows-1 && id_column != 0 && id_column != columns-1)
	{
		for (int i = 1; i < imax - 1; i++)
			for (int j = 1; j < jmax - 1; j++)
				new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];
	}

	// 如果为上边界进程
	if (id_row == 0)
	{
		if (id_column == 0)
		{
			// if (id==0)
			// {
			// 	cout << "p0迭代" << endl;
			// }
			
			for (int i = 2; i < imax - 1; i++)	//不迭代上边界
				for (int j = 2; j < jmax - 1; j++)	//不迭代左边界
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else if (id_column == columns-1)
		{
			for (int i = 2; i < imax - 1; i++)	//不迭代上边界
				for (int j = 1; j < jmax - 2; j++)	//不迭代右边界
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else
		{
			for (int i = 2; i < imax - 1; i++)	//不迭代上边界
				for (int j = 1; j < jmax - 1; j++)	//左右边界迭代
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];
		}

	}

	// 如果为下边界进程
	if (id_row == rows-1)
	{
		if (id_column == 0)
		{
			for (int i = 1; i < imax - 2; i++)	//不迭代下边界
				for (int j = 2; j < jmax - 1; j++)	//不迭代左边界
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else if (id_column == columns-1)
		{
			for (int i = 1; i < imax - 2; i++)	//不迭代下边界
				for (int j = 1; j < jmax - 2; j++)	//不迭代右边界
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];

		} else
		{
			for (int i = 1; i < imax - 2; i++)	//不迭代下边界
				for (int j = 1; j < jmax - 1; j++) // 左右边界迭代
					new_grid[i][j] = pow(dt * c, 2.0) * ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) + (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) + 2.0 * grid[i][j] - old_grid[i][j];
		}

	}
	
		// no need to set Neumann boundaries

	// 要根据位置设置边界值
	// 如果是上边界进程
	if(id_row == 0)
	{
		//需要设置上边界Neumann boundaries
		for (int j = 1; j < jmax-1; j++)
		{
			new_grid[1][j] = new_grid[2][j];
		}

	}


	// 如果是下边界进程
	if(id_row == rows-1)
	{
		//需要设置下边界Neumann boundaries
		for (int j = 1; j < jmax-1; j++)
		{
			new_grid[imax-2][j] = new_grid[imax-3][j];
		}
	}

	// 如果是左边界进程
	if(id_column == 0)
	{
			//设置boundry值（Neumann boundaries）
		for (int i = 1; i < imax-1; i++)
		{
			new_grid[i][1] = new_grid[i][2];
		}

	}

	// 如果是右边界进程
	if(id_column == columns-1)
	{
		//设置boundry值（Neumann boundaries）
		for (int i = 1; i < imax-1; i++)
		{
			new_grid[i][jmax - 2] = new_grid[i][jmax - 3];
		}
	}

	t += dt;

	old_grid.swap(new_grid);	//swap函数交换的是指针，更快
	old_grid.swap(grid);
}


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);	//进程id
	MPI_Comm_size(MPI_COMM_WORLD, &p);	//进程数

	find_dimensions(p, rows, columns);	//根据p大小生成进程结构, 赋值rows和columns的值 rows, columns

	id_to_index(id, id_row, id_column);	//通过进程的id得到其grid行列index， 赋值该进程的进程结构中的行列索引 id_row, id_column

	int grid_rows, grid_columns; //平均(包括ghost layer)

	//分配每个进程的domain
	allocate_grid(grid, new_grid, old_grid, imax, jmax, grid_rows, grid_columns);

	dx = x_max / ((double)imax - 1);	//delta x
	dy = y_max / ((double)jmax - 1);	//delta y

	t = 0.0;

	dt = 0.1 * min(dx, dy) / c; //由于这是一个显式解，因此有严格的时间步长约束， delta t

	int out_cnt = 0, it = 0;	//输出图像的index，迭代的index

	init(grid, old_grid, dx, dy, id_row, id_column, grid_rows, grid_columns);

	grid_to_file(id, out_cnt, grid);
	out_cnt++;
	t_out += dt_out;

	//发送数据缓存
	double** data_array;
	double** recv_array;


	// 构建发送棋盘数组和接受数组内存
	int rows_remain = imax % rows;
	int columns_remain = jmax % columns;

    //使用该进程grid创建发送棋盘
	const int m = grid.size()-2;
	const int n = grid[0].size()-2;
	data_array = new double*[m];
	for (int i = 0; i < m; i++)
		data_array[i] = new double[n];

    //这里接受棋盘按最大的gird来建立，防止截断误差
	const int l = grid_rows+rows_remain;
	const int k = grid_columns+columns_remain;
	recv_array = new double*[l];
	for (int i = 0; i < l; i++)
		recv_array[i] = new double[k];


	// 创建自定义数据类型
	createdatatypes(data_array, m, n);
	createdatatypes_recv(recv_array, l, k);

	while (t < t_max)
	{
		
		//准备发送数据
        int gr = grid.size() - 2; // 该grid有效行数(非ghost)
        int gc = grid[0].size() - 2;    //该grid有效列数(非ghost)
		for (int i = 0; i < gr; i++)
			for (int j = 0; j < gc; j++)
				data_array[i][j] =(double) grid[i+1][j+1]; //将grid数据复制到buffer上，buffer有效可存数据size比gird有效数据大

		for (int i = 0; i < l; i++)
			for (int j = 0; j < k; j++)
				recv_array[i][j] = -1;

		//发送
		int cnt = 0;	//邻居发送计数

		MPI_Request* request = new MPI_Request[4 * 2];	//存储请求:上左右下

		bool recv_top = false;
		bool recv_left = false;
		bool recv_right = false;
		bool recv_bottom = false;
		
		// 对周围进程进行迭代
		for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
			{
				//要发送的grid的行列index（即要发送信息给此grid）
				int com_i = id_row + i;
				int com_j = id_column + j;
				// cout << "p:" << id << "com_i: " << com_i << "com_j: " << com_j;

				//要发送的grid的所属进程id
				int com_id = id_from_index(com_i, com_j);

				//如果该grid的id不为当前进程id（即不属于当前进程的grid），且该grid是合法grid（即属于所有grid里的一个）
				if (com_id != id && com_id >= 0 && com_id < p)
				{
					//上边进程
					if (com_i == id_row-1 && com_j == id_column)
					{
						// tag_num为迭代index
						// cout << "发送给上" << endl;
						MPI_Isend(data_array, 1, Datatype_top, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_top_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						recv_top = true;
						cnt++;
					} 
					else if (com_i == id_row && com_j == id_column-1) //左边进程
					{
						// cout << "发送给左" << endl;
						MPI_Isend(data_array, 1, Datatype_left, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_left_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						recv_left = true;
						cnt++;
					} 
					else if (com_i == id_row && com_j == id_column+1)	//右边进程
					{
						// cout << "发送给右" << endl;
						MPI_Isend(data_array, 1, Datatype_right, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_right_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						recv_right = true;
						cnt++;
					}
					else if (com_i == id_row+1 && com_j == id_column) //下边进程
					{
						// cout << "发送给下" << endl;
						MPI_Isend(data_array, 1, Datatype_bottom, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
						MPI_Irecv(recv_array, 1, Datatype_bottom_recv, com_id, it, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
						recv_bottom = true;
						cnt++;
					}

				}
			}

		//接受
		// cout << "P" << id << "发送邻居: " << cnt << endl;
		MPI_Waitall(cnt * 2, request, MPI_STATUS_IGNORE);

		//将边界值赋值
		copy_boundry(grid, l, k, recv_array, recv_left, recv_top, recv_right, recv_bottom);

		// 计算

		do_iteration();

		if (t_out <= t)
		{
			cout << "P" << id << "\toutput: " << out_cnt << "\tt: " << t << "\titeration: " << it << endl;
			grid_to_file(id, out_cnt, grid);
			out_cnt++;
			t_out += dt_out;	//dt_out = 0.04 打印一次
		}

		it++;

		delete[] request;
	}


	MPI_Finalize();

	delete[] data_array;
	delete[] recv_array;
}
