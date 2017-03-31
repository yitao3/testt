// testt.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "iostream"
#include "testt.h"
#include "cstring"
#include "fstream"
#include "string"
#include <Eigen/Sparse>  
#include <Eigen/Cholesky>
#include <windows.h>
#include<time.h>
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using namespace std;

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;


const int ROW_t = 5687;
const int VOL_t = 3;
const int ROW_c = 2945;
const int VOL_c = 3;
const int ROW_r = 2945;
const int VOL_r = 9;
const int ROW_m = 2945;
const int VOL_m = 8;
/*
const int ROW_t = 223;
const int VOL_t = 3;
const int ROW_c = 133;
const int VOL_c = 3;
const int ROW_r = 133;
const int VOL_r = 9;
const int ROW_m = 133;
const int VOL_m = 8;
*/
int _tmain(int argc, _TCHAR* argv[])
{
	VectorXd cal_coefficient(int node_index, int **relation_mid, element *p, int **relation, double **coord);
	VectorXi fill_node(int start_node, int **relation, int **relation_mid, element *p, double **coord, VectorXd &fill_time);
	
	string  to_s="topo_data_a.txt";// 单元含有哪些节点
	string  co_s = "node_data_a.txt";//[2.33 3.22 1]最后一位判断边界
	string re_s = "relation_data_a.txt";//[2 3 4 0 0 0 1]最后一位是判断边界
	string  mid_s = "relation_mid_data_a.txt";//[3, 11 ,14,0,0,0]存在于哪些单元中
	/*
	string  to_s = "topo_data.txt";// 单元含有哪些节点
	string  co_s = "node_data.txt";//[2.33 3.22 1]最后一位判断边界
	string re_s = "relation_data.txt";//[2 3 4 0 0 0 1]最后一位是判断边界
	string  mid_s = "relation_mid_data.txt";//[3, 11 ,14,0,0,0]存在于哪些单元中
	*/
	int **topo;
	double **coord;
	int **relation;
	int **relation_mid;

	topo = new int*[ROW_t];    //注意，int*[10]表示一个有10个元素的指针数组
	for (int i = 0; i != ROW_t; ++i)
	{
		topo[i] = new int[VOL_t];
	}

	coord = new double*[ROW_c];    //注意，int*[10]表示一个有10个元素的指针数组
	for (int i = 0; i != ROW_c; ++i)
	{
		coord[i] = new double[VOL_c];
	}

	relation = new int*[ROW_r];    //注意，int*[10]表示一个有10个元素的指针数组
	for (int i = 0; i != ROW_r; ++i)
	{
		relation[i] = new int[VOL_r];
	}

	relation_mid = new int*[ROW_m];    //注意，int*[10]表示一个有10个元素的指针数组
	
	for (int i = 0; i != ROW_m; ++i)
	{
		relation_mid[i] = new int[VOL_m];
	}
	topo = read_data(to_s, ROW_t, VOL_t);
	coord = read_data_d(co_s, ROW_c, VOL_c);
	relation = read_data(re_s, ROW_r, VOL_r-1);
	relation_mid = read_data(mid_s, ROW_m, VOL_m-1);

	for (int i = 0; i <= ROW_r-1; i++){
		int k1 = 0;
		int k2 = 0;
		for (int j = 0; j <= VOL_r - 3; j++)
		{

			if (relation[i][j] != 0){
			k1++;
			}
			if (relation_mid[i][j] != 0){
				k2++;
			}
		}
		relation[i][VOL_r - 1] = k1;
		relation_mid[i][VOL_m - 1] = k2;

	}
	element *p;
	p = getelement(coord, topo, relation_mid, ROW_c, ROW_t);

	/*
	for (int i = 0; i <= 30; i++)
	{
		for (int j = 0; j <= 2; j++)
		{
			cout << (p + i)->node_index[j] << "	" << (p + i)->node_mid[j] << "	" << (p + i)->tri_mid[j] << endl;
		}
		cout << endl;
	}
	int i;
	cin >> i;

	MatrixXd m(2, 2);
	m(0, 0) = 1;
	m(1, 0) = 2;
	m(0, 1) = 3;
	m(1, 1) = 4;
	VectorXd b(3);
	b(0) = 1;
	b(1) = 2;
	b(2) = 3;
    VectorXd c = b.transpose();
	VectorXd d;

	cout << -b << endl << b;
	cout << d<<endl<<m*m<<endl<<b.dot(b);
	int hehe;
	cin >> hehe;
	*/
	Vector3i cal_node(32, 132, 22);
	int active = 1;

	VectorXd result(cal_node.size());
	//cout << cal_node << endl << endl;
	//cout << element_index << endl << endl;
	//cout << node_adjacent << endl << endl;
	//VectorXd cal_coefficient(int node_index, int **relation_mid, element *p, int **relation, double **coord);
	VectorXi fill_sequence = VectorXi::Zero(ROW_c);
	VectorXd fill_time = VectorXd::Zero(ROW_c);
	fill_time(0) = 0.05;//初始节点填充时间
	fill_sequence=fill_node(123, relation, relation_mid, p ,coord, fill_time);
	ofstream fout("output.txt");
	for (int i = 0; i < ROW_c; i++)
	{
		fout << fill_sequence(i)<<"	"<<fill_time(i)<<endl;
	}
	return 0;
}
VectorXd cal_coefficient(int node_index, int **relation_mid, element *p, int **relation, double **coord)
{
	if (node_index == 15)
	{
		int a = 10;
	}
	int k = relation_mid[node_index][VOL_m - 1];
	int k2 = relation[node_index][VOL_r - 1];
	VectorXi element_index(k);
	VectorXi node_adjacent(k2);
	for (int i = 0; i <= k - 1; i++){
		int index = node_index;
		element_index(i) = relation_mid[index][i];
	}
	for (int i = 0; i <= k2 - 1; i++){
		int index = node_index;
		node_adjacent(i) = relation[index][i];
	}

	int ele_num = element_index.size();//警惕：有时候周围节点数和周围单元数不同
	VectorXd co = VectorXd::Zero(k2 + 1);
	VectorXi sequence(k2 + 1);//包括主节点和周围节点的序列
	sequence(0) = node_index;

	for (int i = 1; i <= k2; i++)
	{
		sequence(i) = node_adjacent(i - 1);
	}
	for (int i = 0; i <= ele_num - 1; i++) //该节点所在的element的每个进行循环
	{
		int ele_i = element_index(i);//该节点所在的element的index
		element *temp = (p + ele_i);
		double area = temp->area;
		double x1 = temp->x1;
		double x2 = temp->x2;
		double x3 = temp->x3;
		double y1 = temp->y1;
		double y2 = temp->y2;
		double y3 = temp->y3;
		Vector3d delta_px(1.0 / 2.0 / area*(y2 - y3), 1.0 / 2.0 / area*(y3 - y1), 1.0 / 2.0 / area*(y1 - y2));
		Vector3d delta_py(1.0 / 2.0 / area*(x3 - x2), 1.0 / 2.0 / area*(x1 - x3), 1.0 / 2.0 / area*(x2 - x1));
		double permeability_x = 1.0;
		double permeability_y = 1.0;
		double viscosity = 1.0;
		int index_123 = 0;
		for (int j = 0; j <= 2; j++)
		{
			if (node_index == temp->node_index[j])
			{
				index_123 = j;

			}
		}
		Vector3d q_total(0, 0, 0);
		Vector3d qx(0, 0, 0);
		Vector3d qy(0, 0, 0);
		for (int j = 0; j <= 2; j++)
		{
			Vector2d mid2_active_node(coord[temp->node_index[index_123]][0] - temp->tri_mid[0], coord[temp->node_index[index_123]][1] - temp->tri_mid[1]);
			if (j != index_123)
			{
				Vector2d mid2mid(temp->node_mid[j] - temp->tri_mid[0], temp->node_mid[j + 3] - temp->tri_mid[1]);//中点到中点
				Vector2d vertical(mid2mid(1), -mid2mid(0));
				if (vertical.dot(mid2_active_node) > 0)
				{
					vertical = -vertical;
				}

				qx = qx + delta_px*vertical(0);
				qy = qy + delta_py*vertical(1);
			}
		}
		q_total = qx*permeability_x / viscosity + qy*permeability_y / viscosity;
		for (int m = 0; m < k2 + 1; m++)
		{
			for (int n = 0; n <= 2; n++)
			{
				if (temp->node_index[n] == sequence(m)){
					co(m) = co(m) + q_total(n);
				}
			}
		}
	}
	//cout << co << endl;
	return -co;
}
VectorXi fill_node(int start_node, int **relation, int **relation_mid, element *p, double **coord, VectorXd &fill_time)
{
	VectorXd cal_pressure(SpMat coefficients, int numFrontnode, VectorXd boundary);
	void cal_q(MatrixXd &flow_rate, VectorXi fill_sequence, VectorXd pressure, VectorXi front_node, int **relation_mid, element *p, int **relation, int numFilled, int numFront, double **coord);
	int node_num = ROW_c;
	VectorXi fill_sequence = VectorXi::Zero(node_num);//完全充填的顺序
	VectorXd pressure = VectorXd::Zero(node_num);//
	MatrixXd flow_rate = MatrixXd::Zero(node_num,5); //流率，里面有的体积以及总的体积,是否被填满,如果被填满那么序列是多少。
	VectorXi front_node = VectorXi::Zero(node_num);//判断哪些是前沿节点
	int numFilled = 1;
	int numFront = 0;
	int num_1 = relation[start_node][VOL_r - 1];
	flow_rate(start_node, 3) = 1;
	flow_rate(start_node, 4) = 1;//这里从1开始
	for (int i = 0; i <=num_1-1 ; i++)//把初始节点的周围节点置于流动前沿vector中
	{
		front_node(relation[start_node][i]) = 1; 
		numFront++;
	}
	for (int i = 0; i <= node_num - 1; i++)//计算每个控制体的体积
	{
		int num_k = relation_mid[i][VOL_m - 1];
		for (int j = 0; j <= num_k - 1; j++)
		{
			flow_rate(i, 2) = flow_rate(i, 2) + 1 / 3.0*((p + i)->area);
		}
			
	}

	fill_sequence(0) = start_node;
	double totalTime = 0;
	double invTime = 0;
	for (int i = 2; i <= node_num-1; i++)//填充进行n-2次
	{


		clock_t start_time_t = clock();
		
		cout << double(i) / node_num<<endl;
		SpMat coefficientMatrix(numFilled, numFilled);
		std::vector < Triplet < double > > triplets;
		clock_t start_time_s = clock();
		
		for (int j = 0; j <= numFilled - 1; j++) //对每个已经填充的节点进行循环
		{
			VectorXd tempCoe = cal_coefficient(fill_sequence[j] , relation_mid,  p , relation, coord);
			//coefficientMatrix(j, j) = tempCoe(0);
			triplets.emplace_back(j, j, tempCoe(0));
		    //cout << tempCoe << endl;
			for (int k = 0; k <= relation[fill_sequence[j]][VOL_m] - 1; k++) //对该填充节点的周围节点进行循环
			{
				int index = relation[fill_sequence[j]][k];
				if (flow_rate(index, 4) != 0)
				{
					triplets.emplace_back(j, int(flow_rate(index, 4))-1, tempCoe(k + 1));  //这里把1减掉
				}
			}
		}

		clock_t end_time_s = clock();
		coefficientMatrix.setFromTriplets(triplets.begin(), triplets.end());    // 初始化系数矩阵
		/*
				for (int j = 0; j <= numFilled - 1; j++) //对每个已经填充的节点进行循环
		{
			VectorXd tempCoe = cal_coefficient(fill_sequence[j] , relation_mid,  p , relation, coord);
			coefficientMatrix(j, j) = tempCoe(0);
		    //cout << tempCoe << endl;
			for (int k = 0; k <= relation[fill_sequence[j]][VOL_m] - 1; k++) //对该填充节点的周围节点进行循环
			{
				for (int l = 0; l <= numFilled - 1; l++)//对该周围节点判断是不是填充满了，是的话就加入参数矩阵中
				{
					if (j == 2)
					{
						int b = 0;
					}
					//cout << "周围的点：" << relation[fill_sequence[j]][k] << endl;
					//cout <<"填满的点："<< fill_sequence(l)<<endl;
					if (relation[fill_sequence[j]][k] == fill_sequence(l))
					{
						if (l == 6)
						{
							int c = 0;
						}
						coefficientMatrix(j, l) = tempCoe(k+1);  //注意要加一
					}
				}
			}
		}
		*/
		//cout << coefficientMatrix << endl;
		VectorXd boundary = VectorXd::Zero(numFilled);
		boundary(0) = 1;//恒定流量 
		VectorXd pressureTemp = VectorXd::Zero(numFilled);

		pressureTemp = cal_pressure(coefficientMatrix, numFilled, boundary);
		//cout << pressureTemp<<endl;
		

		for (int j = 0; j <= numFilled - 1; j++)
		{
			pressure(fill_sequence(j)) = pressureTemp(j);
		}
		cal_q(flow_rate, fill_sequence, pressure, front_node, relation_mid, p, relation, numFilled, numFront,coord);
		double min_time = 100000;
		int min_index = 0;
		numFront = 0;
		for (int j = 1; j <= node_num - 1; j++)
		{
			double c = flow_rate(j, 3);
			double d = flow_rate(j, 0);
			if (flow_rate(j, 3) != 1 && flow_rate(j, 0) != 0)
			{
				//cout << flow_rate(j, 0) << endl;
				//cout <<"前沿节点"<< j << endl;
				front_node(j) = 1;
				numFront++;
				double fillTime = (flow_rate(j, 2) - flow_rate(j, 1)) / flow_rate(j, 0);
				if ( fillTime < min_time)
				{
					min_time = fillTime;
					min_index = j;
				}
			}
		}
		for (int j = 1; j <= node_num - 1; j++)
		{
			if (flow_rate(j, 3) != 1 && flow_rate(j, 0) != 0)
			{
				flow_rate(j, 1) = flow_rate(j, 1)+flow_rate(j, 0)*min_time;
			}
		}
		cout << "新节点:" << min_index << endl;
		if (min_index == 62)
		{
			int d = 0;
		}
		for (int j = 0; j <= relation[min_index][VOL_r - 1] - 1; j++)//把新节点点的周围节点置于流动前沿vector中
		{
			if (flow_rate(relation[min_index][j], 3) != 1 && front_node(relation[min_index][j])!=1) //如果该节点没被填满
			{
				numFront++;
				//cout << "前沿节点" << relation[min_index][j] << endl;
				front_node(relation[min_index][j]) = 1;
			}
			
		}
		fill_sequence(numFilled) = min_index;  //把新节点加入到已充满中
		flow_rate(min_index, 4) = numFilled+1;
		flow_rate(min_index, 3) = 1; //flow_rate中把新节点设为已充满
		front_node(min_index) = 0;   //把新节点从前沿节点中删除
		numFront--;
		/*
		for (int j = 1; j <= node_num - 1; j++)
		{
			if (front_node(j)==1)
			{

				cout << "前沿节点" << j << endl;
			}
		}
		*/
		numFilled++;
		fill_time(numFilled-1) = min_time;
		clock_t end_time_t = clock();
		totalTime += static_cast<double>(end_time_t - start_time_t) / CLOCKS_PER_SEC * 1000;
		invTime += static_cast<double>(end_time_s - start_time_s) / CLOCKS_PER_SEC * 1000;
		cout << "总时间: " << totalTime << endl;
		cout << "分时间: " << invTime << endl;

		//cout << "矩阵求逆时间: " << all_s << endl;

	}
	return fill_sequence;
}
/*
VectorXd cal_pressure(vector<T> &coefficients, int numFrontnode, VectorXd boundary)
{
	SpMat A(numFrontnode, numFrontnode);
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	Eigen::SimplicialCholesky<SpMat> chol(A);
	Eigen::VectorXd x = chol.solve(boundary);
	return x;
}
*/
VectorXd cal_pressure(SpMat coefficientMatrix, int numFrontnode, VectorXd boundary)
{
	Eigen::SimplicialCholesky<SpMat> chol(coefficientMatrix);  // performs a Cholesky factorization of A
	Eigen::VectorXd x = chol.solve(boundary);
	return x;
	//VectorXd x = VectorXd::Zero(numFrontnode);
	//return x = coefficients.ldlt().solve(boundary);
}
void cal_q(MatrixXd &flow_rate, VectorXi fill_sequence, VectorXd pressure, VectorXi front_node, int **relation_mid, element *p, int **relation, int numFilled, int numFront , double **coord)
{
	double tempQ = 0;
	int j = 0;//判断是不是前沿节点
	for (int i = 0; i <= numFront - 1; i++)
	{
		while (front_node[j] == 0)//直到找到下一个前沿节点
		{
			j = j + 1;
		}
		VectorXd coefficientJ = cal_coefficient(j, relation_mid, p, relation, coord);
		int adjacentNum = relation[j][VOL_r-1];
		for (int k = 0; k <= adjacentNum-1; k++)
		{
			tempQ = tempQ + coefficientJ(k+1)*pressure(relation[j][k]);


		}
		
		flow_rate(j, 0) = -tempQ;
		tempQ = 0;
		j = j + 1;
	}

}

//VectorXd cal_coefficient(VectorXi cal_node, int active, VectorXi element_index, element *ele, VectorXi node_adjacent) //element_index 包括节点的所有单元,node_adjacent相邻节点


