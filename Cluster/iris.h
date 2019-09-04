#pragma once
#include"mysql.h"
#include<iostream>
#include<vector>
#include<stack>
#include<cmath>
#include<string>
using namespace std;
#include"dbconnect.h"
class irisnode {
public :
	int id;
	float sepal_length;//萼片
	float sepal_width;
	float petal_length;//花瓣
	float petal_width;
	string category;
	irisnode(int i,float sepal_l, float sepal_w, float petal_l, float petal_w, string cate) {
		id = i;
		sepal_length = sepal_l;
		sepal_width = sepal_w;
		petal_length = petal_l;
		petal_width = petal_w;
		category = cate;
	}
};
float selMax = INT_MIN;
float selMin = INT_MAX;
float sewMax = INT_MIN;
float sewMin = INT_MAX;
float pelMax = INT_MIN;
float pelMin = INT_MAX;
float pewMax = INT_MIN;
float pewMin = INT_MAX;
float CSv= INT_MAX;
class iris
{
public:
	int id;
	float sepal_length;//萼片
	float sepal_width;
	float petal_length;//花瓣
	float petal_width;
	vector<irisnode> nds;
	string category;
	iris();
	iris(int i,float sepal_l, float sepal_w, float petal_l, float petal_w, string cate);
	~iris();
	vector<iris> init_iris();
	float irisDistance(iris i1, iris i2);
	void iris_cluster();
	void iris_cluster2();
	void FACA_DTRS();
private:
	
};
iris::iris(int i,float sepal_l, float sepal_w, float petal_l, float petal_w, string cate) {
	id = i;
	sepal_length = sepal_l;
	sepal_width = sepal_w;
	petal_length = petal_l;
	petal_width = petal_w;
	category = cate;
	irisnode t(id,sepal_length, sepal_width, petal_length, petal_width, category);
	nds.push_back(t);
}
iris::iris() {

}
iris::~iris() {
	
}

float iris::irisDistance(iris i1, iris i2) {
	return (abs(i1.sepal_length - i2.sepal_length) / (selMax - selMin) +
		abs(i1.sepal_width-i2.sepal_width)/(sewMax-sewMin)+
		abs(i1.petal_length-i2.petal_length)/(pelMax-pelMin)+
		abs(i1.petal_width-i2.petal_width)/(pewMax-pewMin)) / 4;

}
//从数据库中拿数据
vector<iris> iris::init_iris() {
	dbconnect *dbcon=new dbconnect();
	MYSQL * con=dbcon->get_con();
	MYSQL_RES *res;
	MYSQL_ROW row;
	char tmp[400];
	char tablename[50] = "iris";
	char * query = NULL;
	int rt;//return value;

	vector<iris> ary;

	sprintf_s(tmp, "select * from %s where ID>=1 and ID<=150", tablename);
	rt = mysql_real_query(con, tmp, strlen(tmp));
	if (rt)
	{
		printf("Error making query: %s !!!\n", mysql_error(con));
	}
	else
	{
		//printf("%s executed!!!\n", tmp);
	}
	res = mysql_store_result(con);//将结果保存在res结构体中

	while (row = mysql_fetch_row(res)) {
		iris crount(atof(row[0]),atof(row[1]), atof(row[2]), atof(row[3]), atof(row[4]),row[5]);
		if (selMax < atof(row[1]))
			selMax = atof(row[1]);
		if (selMin > atof(row[1]))
			selMin = atof(row[1]);
		if (sewMax < atof(row[2]))
			sewMax = atof(row[2]);
		if (sewMin > atof(row[2]))
			sewMin = atof(row[2]);
		if (pelMax < atof(row[3]))
			pelMax = atof(row[3]);
		if (pelMin > atof(row[3]))
			pelMin = atof(row[3]);
		if (pewMax < atof(row[4]))
			pewMax = atof(row[4]);
		if (pewMin > atof(row[4]))
			pewMin = atof(row[4]);
		ary.push_back(crount);
	}
	mysql_free_result(res);
	free(dbcon);
	return ary;
}

void iris::iris_cluster() {
	vector<iris> irisSet;
	vector<vector<float>> dTable;
	irisSet = init_iris();
	for (int i = 0; i < irisSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < irisSet.size(); j++) {
			if (j > i)
				temp.push_back(irisDistance(irisSet[i], irisSet[j]));
			else if (j < i)
				temp.push_back(dTable[j][i]);
			else
				temp.push_back(0);
		}
		dTable.push_back(temp);
	}

	while (true) {
		float mindt = INT_MAX;
		int mi, mj,k=0,flag=0;
		float val=0,CSn=0;

		for (int i = 0; i < dTable.size(); i++) {
			for (int j = i+1; j < dTable[i].size(); j++) {
				val += dTable[i][j]; k++;
				if (dTable[i][j] != dTable[0][1])
					flag = 1;
				if (dTable[i][j] < mindt) {
					mindt = dTable[i][j];
					mi = i; mj = j;
				}
			}
		}

		if (flag == 0) {
			cout << irisSet.size() << endl;
			return;
		}

		val =1- val  / k;

		for (int i = 0; i < dTable.size(); i++) {
			for (int j = i+1; j < dTable[i].size(); j++) {
				if (1 - dTable[i][j] >= val)
					CSn += (0.5 - ((1 - dTable[i][j]) - val) / (2 - 2 * val));
				if(1-dTable[i][j]<val)
					CSn += (0.5 + (val-(1 - dTable[i][j]) ) / (2 * val));
			}
		}
		
		if (CSv > CSn&&irisSet.size() == 1) {
			cout << 1 << endl;
			return;
		}

		if (CSv <= CSn) {
			cout << irisSet.size() + 1<<endl;
			return;
		}

		CSv = CSn;

		for (int i = 0; i < irisSet[mj].nds.size(); i++) {
			irisSet[mi].nds.push_back(irisSet[mj].nds[i]);
		}

		vector<iris>::iterator imj = irisSet.begin();
		imj += mj;
		irisSet.erase(imj);
		//最近邻聚类
		for (int j = 0; j < dTable.size(); j++) {
			if (j == mi || j == mj) continue;
			if (dTable[mi][j] > dTable[mj][j]) {
				dTable[mi][j] = dTable[mj][j];
				dTable[j][mi] = dTable[mi][j];
			}
		}

		vector<vector<float>>::iterator ii = dTable.begin();
		ii += mj;
		dTable.erase(ii);
		for (int i = 0; i < dTable.size(); i++) {
			vector<float>::iterator ij = dTable[i].begin();
			ij += mj;
			dTable[i].erase(ij);
		}
	}
	/*for (int i = 0; i < irisSet.size(); i++) {
		vector<node>::iterator it = irisSet[i].nds.begin();
		for (int j = 0; j < irisSet[i].nds.size(); j++) {
			cout << (*it).sepal_length << "  ";
			cout << (*it).sepal_width << "  ";
			cout << (*it).petal_length<<"  ";
			cout << (*it).petal_width<<"  ";
			cout<<(*it).category<<endl;
			it++;
		}
		cout << "-------------------";
		cout << endl;
	}
	cout << irisSet[0].nds.size() << "  ";
	cout << irisSet[1].nds.size() << "  ";*/
}

void iris::iris_cluster2() {
	vector<iris> irisSet;
	vector<vector<float>> dTable;
	irisSet = init_iris();
	for (int i = 0; i < irisSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < irisSet.size(); j++) {
			if (j > i)
				temp.push_back(irisDistance(irisSet[i], irisSet[j]));
			else if (j < i)
				temp.push_back(dTable[j][i]);
			else
				temp.push_back(0);
		}
		dTable.push_back(temp);
	}

	float val = 0;
	int k = 0;

	for (int i = 0; i < dTable.size(); i++) {
		for (int j = i + 1; j < dTable[i].size(); j++) {
			val += dTable[i][j]; k++;
		}
	}

	val = 1-val / k;

	vector<vector<float>> P;

	for (int i = 0; i < dTable.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < dTable[i].size(); j++) {
			if (j > i) {
				if ((1 - dTable[i][j]) >= val)
					temp.push_back(0.5 + ((1 - dTable[i][j] - val) / (2 - 2 * val)));
				else
					temp.push_back(0.5 - (val - (1 - dTable[i][j])) / (2 * val));
			}
				
			else if (j < i)
				temp.push_back(P[j][i]);
			else
				temp.push_back(0);
		}
		P.push_back(temp);
	}

	while (true) {
		float maxf = INT_MIN;
		int mi=0, mj=0;
		vector<vector<float>> f;
		f.clear();

		for (int i = 0; i < irisSet.size(); i++) {
			vector<float> temp;
			for (int j = 0; j <irisSet.size(); j++) {
				if (j > i) {
					float sump = 0;
					for (vector<irisnode>::iterator it1 = irisSet[i].nds.begin(); it1 < irisSet[i].nds.end();it1++)
						for(vector<irisnode>::iterator it2 = irisSet[j].nds.begin(); it2 < irisSet[j].nds.end(); it2++)
							sump+=P[(*it1).id-1][(*it2).id-1];
					temp.push_back(sump / (irisSet[i].nds.size()*irisSet[j].nds.size()));
				}
				else if (j < i)
					temp.push_back(f[j][i]);
				else
					temp.push_back(0);
			}
			f.push_back(temp);
		}

		for (int i = 0; i < f.size(); i++) {
			for (int j = i + 1; j < f[i].size(); j++) {
				if (f[i][j] > maxf) {
					maxf = f[i][j];
					mi = i; mj = j;
				}
			}
		}

		if (maxf <= 0.5) {
			cout << f.size()<<endl;
			return;
		}

		if (maxf > 0.5&&irisSet.size() == 1) {
			cout << "error!" << endl;
			return;
		}

		for (vector<irisnode>::iterator it1 = irisSet[mi].nds.begin(); it1 < irisSet[mi].nds.end(); it1++)
			for (vector<irisnode>::iterator it2 = irisSet[mj].nds.begin(); it2 < irisSet[mj].nds.end(); it2++) {
				P[(*it1).id-1][(*it2).id-1] = 1 - P[(*it1).id-1][(*it2).id-1];
			}

		for (int i = 0; i < irisSet[mj].nds.size(); i++) {
			irisSet[mi].nds.push_back(irisSet[mj].nds[i]);
		}

		vector<iris>::iterator imj = irisSet.begin();
		imj += mj;
		irisSet.erase(imj);

	}
}

void iris::FACA_DTRS() {
	/*vector<iris> irisSet;
	vector<vector<float>> dTable;
	irisSet = init_iris();
	for (int i = 0; i < irisSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < irisSet.size(); j++) {
			if (j > i)
				temp.push_back(irisDistance(irisSet[i], irisSet[j]));
			else if (j < i)
				temp.push_back(dTable[j][i]);
			else
				temp.push_back(0);
		}
		dTable.push_back(temp);
	}

	float val = 0;
	int k = 0;

	for (int i = 0; i < dTable.size(); i++) {
		for (int j = i + 1; j < dTable[i].size(); j++) {
			val += dTable[i][j]; k++;
		}
	}

	val = 1 - val / k;

	vector<vector<float>> P;

	for (int i = 0; i < dTable.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < dTable[i].size(); j++) {
			if (j > i) {
				if ((1 - dTable[i][j]) >= val)
					temp.push_back(0.5 + ((1 - dTable[i][j] - val) / (2 - 2 * val)));
				else
					temp.push_back(0.5 - (val - (1 - dTable[i][j])) / (2 * val));
			}

			else if (j < i)
				temp.push_back(P[j][i]);
			else
				temp.push_back(0);
		}
		P.push_back(temp);
	}

	while (true) {
		float maxf = INT_MIN;
		vector<vector<float>> f;
		f.clear();

		for (int i = 0; i < irisSet.size(); i++) {
			vector<float> temp;
			for (int j = 0; j < irisSet.size(); j++) {
				if (j > i) {
					float sump = 0;
					for (vector<node>::iterator it1 = irisSet[i].nds.begin(); it1 < irisSet[i].nds.end(); it1++)
						for (vector<node>::iterator it2 = irisSet[j].nds.begin(); it2 < irisSet[j].nds.end(); it2++)
							sump += P[(*it1).id - 1][(*it2).id - 1];
					temp.push_back(sump / (irisSet[i].nds.size()*irisSet[j].nds.size()));
				}
				else if (j < i)
					temp.push_back(f[j][i]);
				else
					temp.push_back(0);
			}
			f.push_back(temp);
		}

		int k2 = sqrt(f.size()),k1=f.size();
		stack<int> mis, mjs;
		while (!mis.empty())mis.pop();
		while (!mjs.empty())mjs.pop();

		for (int i = 0; i < f.size(); i++) {
			for (int j = i + 1; j < f[i].size(); j++) {
				if (f[i][j] > maxf) {
					maxf = f[i][j];
					mis.push(i);
					mjs.push(j);
				}
			}


		}

		while (!mis.empty()) {
			int mi = mis.top();
			int mj = mjs.top();
			mis.pop(); mjs.pop();
			for (vector<node>::iterator it1 = irisSet[mi].nds.begin(); it1 < irisSet[mi].nds.end(); it1++)
				for (vector<node>::iterator it2 = irisSet[mj].nds.begin(); it2 < irisSet[mj].nds.end(); it2++) {
					P[(*it1).id - 1][(*it2).id - 1] = 1 - P[(*it1).id - 1][(*it2).id - 1];
				}

			for (int i = 0; i < irisSet[mj].nds.size(); i++) {
				irisSet[mi].nds.push_back(irisSet[mj].nds[i]);
			}

			vector<iris>::iterator imj = irisSet.begin();
			imj += mj;
			irisSet.erase(imj);
		}
	}*/
}

