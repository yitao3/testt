#include <iostream>
#include "fstream"
using namespace std;

struct element
{
	int node_index[3];
	double node_mid[6];
	double tri_mid[2];
	double area;
	double x1;
	double x2;
	double x3;
	double y1;
	double y2;
	double y3;

};


element *getelement(double **coord, int **topo, int **relation_mid, int node_num, int ele_num){
	double area_tri(double x1, double y1, double x2, double y2, double x3, double y3);
	element *p;
	p = new element[ele_num];
	int i;
	int node_i[3];
	double node_m[6];
	double tri_m[2] = { 0, 0 };
	int size_node_i = sizeof(node_i);
	int size_node_m = sizeof(node_m);
	int size_tri_m = sizeof(tri_m);
	for (i = 0; i < ele_num; i++)
	{
		//cout << topo[i][0] << endl;
		//int a = topo[i][0];
		double x1 = coord[topo[i][0]][0];
		double x2 = coord[topo[i][1]][0];
		double x3 = coord[topo[i][2]][0];
		double y1 = coord[topo[i][0]][1];
		double y2 = coord[topo[i][1]][1];
		double y3 = coord[topo[i][2]][1];
		tri_m[0] = (x1+x2+x3) / 3.0;
		tri_m[1] = (y1+y2+y3) / 3.0;
		node_i[0] = topo[i][0];
		node_i[1] = topo[i][1];
		node_i[2] = topo[i][2];
		int j;
		for (j = 0; j <= 2; j++)
		{
			node_m[j] = (tri_m[0] * 3 - coord[int(topo[i][j])][0]) / 2.0;
			node_m[j + 3] = (tri_m[1] * 3 - coord[int(topo[i][j])][1]) / 2.0;
		}
		memcpy((p + i)->node_index, node_i, size_node_i);
		memcpy((p + i)->node_mid, node_m, size_node_m);
		memcpy((p + i)->tri_mid, tri_m, size_tri_m);
		(p + i)->area = area_tri(x1, y1, x2, y2, x3, y3);
		(p + i)->x1 = x1;
		(p + i)->x2 = x2;
		(p + i)->x3 = x3;
		(p + i)->y1 = y1;
		(p + i)->y2 = y2;
		(p + i)->y3 = y3;
	}


	return p;
}
int **read_data(string f_o, int row, int column)
{
	ifstream f(f_o);
	int **p;
	/*
	if (!f)
	{
		cout << "hehe";
	}
	else{
		cout << "you";
	}
	*/
	p = new int*[row];    //注意，int*[10]表示一个有10个元素的指针数组
	for (int i = 0; i != row; ++i)
	{
		p[i] = new int[column];
	}
	for (int i = 0; i < row; ++i){
		for (int j = 0; j < column; ++j){
			f >> p[i][j];
		}
	}
	f.close();
	return p;
}
double area_tri(double x1, double y1, double x2, double y2, double x3, double y3)
{
	return 0.5*(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
}

double **read_data_d(string f_o, int row, int column)
{
	ifstream f(f_o);
	double **p;
	p = new double*[row];    //注意，int*[10]表示一个有10个元素的指针数组
	for (int i = 0; i != row; ++i)
	{
		p[i] = new double[column];
	}
	for (int i = 0; i < row; ++i){
		for (int j = 0; j < column; ++j){
			f >> p[i][j];
		}
	}
	f.close();
	return p;
}