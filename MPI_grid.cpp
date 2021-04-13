#include <mpi.h>
#include <iostream>

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


int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);	//进程id
	MPI_Comm_size(MPI_COMM_WORLD, &p);	//进程数


	find_dimensions(p, rows, columns);	//根据p大小生成gird

	id_to_index(id, id_row, id_column);	//通过进程的id得到其grid行列index

	int *data_to_send = new int[8 * 3];
	int *data_to_recv = new int[8 * 3];

	int cnt = 0;	//邻居发送计数

	MPI_Request* request = new MPI_Request[8 * 2];	//存储请求


	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
		{
			//要发送的grid（即要发送信息给此grid）
			int com_i = id_row + i;
			int com_j = id_column + j;

			//要发送的grid的所属进程id
			int com_id = id_from_index(com_i, com_j);

			//如果该grid的id不为当前进程id（即不属于当前进程的grid），且该grid是合法grid（即属于所有grid里的一个）
			if (com_id != id && com_id >= 0 && com_id < p)
			{
				//存储当前进程的id，以及进程分配的grid行列index
				data_to_send[cnt * 3] = id;
				data_to_send[cnt * 3 + 1] = id_row;
				data_to_send[cnt * 3 + 2] = id_column;

				//当前进程发送自己的行列index和id给要发送的gird所属进程
				MPI_Isend(&data_to_send[cnt * 3], 3, MPI_INT, com_id, tag_num, MPI_COMM_WORLD, &request[cnt * 2]);
				MPI_Irecv(&data_to_recv[cnt * 3], 3, MPI_INT, com_id, tag_num, MPI_COMM_WORLD, &request[cnt * 2 + 1]);
				cnt++;
			}
		}

	MPI_Waitall(cnt * 2, request, MPI_STATUS_IGNORE);

	//每个进程输出：当前进程id：收到的所有邻居id，以及对应的邻居index
	for (int i = 0; i < cnt; i++)
	{
		cout << id << " from " << data_to_recv[i * 3] << " ( " << data_to_recv[i * 3 + 1] << ", " << data_to_recv[i * 3 + 2] << ")" << endl;
	}

	MPI_Finalize();

	delete[] data_to_send;
	delete[] data_to_recv;
	delete[] request;
}
