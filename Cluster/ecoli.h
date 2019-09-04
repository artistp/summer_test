#pragma once
#include"mysql.h"
#include<iostream>
#include<vector>
#include<stack>
#include<cmath>
#include<string>
using namespace std;
#include"dbconnect.h"
class ecolinode {
public:
	int id;
	string Sequence_Name;
	float mcg;
	float gvh;
	float lip;
	float chg;
	float aac;
	float alm1;
	float alm2;
	string category;
	ecolinode(int i,string Sequence_Nam,float mc,float gv,float li,float ch,float aa,float al1,float al2,string categor) {
		id=i;
		Sequence_Name=Sequence_Nam;
		mcg=mc;
		gvh=gv;
		lip=li;
		chg=ch;
		aac=aa;
		alm1=al1;
		alm2=al2;
		category=categor;
	}
};
float mcgMax = INT_MIN;
float mcgMin = INT_MAX;
float gvhMax = INT_MIN;
float gvhMin = INT_MAX;
float lipMax = INT_MIN;
float lipMin = INT_MAX;
float chgMax = INT_MIN;
float chgMin = INT_MAX;
float aacMax = INT_MIN;
float aacMin = INT_MAX;
float alm1Max = INT_MIN;
float alm1Min = INT_MAX;
float alm2Max = INT_MIN;
float alm2Min = INT_MAX;
class ecoli
{
public:
	int id;
	string Sequence_Name;
	float mcg;
	float gvh;
	float lip;
	float chg;
	float aac;
	float alm1;
	float alm2;
	string category;
	vector<ecolinode> nds;
	ecoli();
	~ecoli();
	ecoli(int i, string Sequence_Nam, float mc, float gv, float li, float ch, float aa, float al1, float al2, string categor);
	vector<ecoli> init_ecoli();
	float ecoliDistance(ecoli e1, ecoli e2);
	void ecoli_cluster();
	void ecoli_cluster2();
};

ecoli::ecoli(int i, string Sequence_Nam, float mc, float gv, float li, float ch, float aa, float al1, float al2, string categor) {
	id = i;
	Sequence_Name = Sequence_Nam;
	mcg = mc;
	gvh = gv;
	lip = li;
	chg = ch;
	aac = aa;
	alm1 = al1;
	alm2 = al2;
	category = categor;
	ecolinode t(id ,Sequence_Name ,mcg ,gvh,lip ,chg ,aac ,alm1 ,alm2 ,category);
	nds.push_back(t);
}
ecoli::ecoli() {

}
ecoli::~ecoli() {

}

float ecoli::ecoliDistance(ecoli e1, ecoli e2) {
	return (abs(e1.mcg - e2.mcg) / (mcgMax - mcgMin) +
		abs(e1.gvh - e2.gvh) / (gvhMax - gvhMin) +
		abs(e1.lip - e2.lip) / (lipMax - lipMin) +
		abs(e1.chg - e2.chg) / (chgMax - chgMin)+
		abs(e1.aac - e2.aac) / (aacMax - aacMin)+
		abs(e1.alm1 - e2.alm1) / (alm1Max - alm1Min)+
		abs(e1.alm2 - e2.alm2) / (alm2Max - alm2Min)) / 7;
}

vector<ecoli> ecoli::init_ecoli() {
	dbconnect *dbcon = new dbconnect();
	MYSQL * con = dbcon->get_con();
	MYSQL_RES *res;
	MYSQL_ROW row;
	char tmp[400];
	char tablename[50] = "ecoli";
	char * query = NULL;
	int rt;//return value;

	vector<ecoli> ary;

	sprintf_s(tmp, "select * from %s", tablename);
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
		ecoli crount(atof(row[0]), row[1], atof(row[2]), atof(row[3]), atof(row[4]), atof(row[5]), atof(row[6]), atof(row[7]), atof(row[8]), row[9]);
		if (mcgMax < atof(row[2]))
			mcgMax = atof(row[2]);
		if (mcgMin > atof(row[2]))
			mcgMin = atof(row[2]);
		if (gvhMax < atof(row[3]))
			gvhMax = atof(row[3]);
		if (gvhMin > atof(row[3]))
			gvhMin = atof(row[3]);
		if (lipMax < atof(row[4]))
			lipMax = atof(row[4]);
		if (lipMin > atof(row[4]))
			lipMin = atof(row[4]);
		if (chgMax < atof(row[5]))
			chgMax = atof(row[5]);
		if (chgMin > atof(row[5]))
			chgMin = atof(row[5]);
		if (aacMax < atof(row[6]))
			aacMax = atof(row[6]);
		if (aacMin > atof(row[6]))
			aacMin = atof(row[6]);
		if (alm1Max < atof(row[7]))
			alm1Max = atof(row[7]);
		if (alm1Min > atof(row[7]))
			alm1Min = atof(row[7]);
		if (alm2Max < atof(row[8]))
			alm2Max = atof(row[8]);
		if (alm2Min > atof(row[8]))
			alm2Min = atof(row[8]);
		ary.push_back(crount);
	}
	mysql_free_result(res);
	free(dbcon);
	return ary;
}

void ecoli::ecoli_cluster() {
	vector<ecoli> ecoliSet;
	vector<vector<float>> dTable;
	ecoliSet = init_ecoli();
	for (int i = 0; i < ecoliSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < ecoliSet.size(); j++) {
			if (j > i)
				temp.push_back(ecoliDistance(ecoliSet[i], ecoliSet[j]));
			else if (j < i)
				temp.push_back(dTable[j][i]);
			else
				temp.push_back(0);
		}
		dTable.push_back(temp);
	}

	while (true) {
		float mindt = INT_MAX;
		int mi, mj, k = 0, flag = 0;
		float val = 0, CSn = 0;
		//
cout << dTable.size();
		for (int i = 0; i < dTable.size(); i++) {
			for (int j = i + 1; j < dTable[i].size(); j++) {
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
			cout << ecoliSet.size() << endl;
			return;
		}

		val = 1 - val / k;

		for (int i = 0; i < dTable.size(); i++) {
			for (int j = i + 1; j < dTable[i].size(); j++) {
				if (1 - dTable[i][j] >= val)
					CSn += (0.5 - ((1 - dTable[i][j]) - val) / (2 - 2 * val));
				if (1 - dTable[i][j] < val)
					CSn += (0.5 + (val - (1 - dTable[i][j])) / (2 * val));
			}
		}

		if (CSv > CSn&&ecoliSet.size() == 1) {
			cout << 1 << endl;
			return;
		}

		if (CSv <= CSn) {
			cout << ecoliSet.size() + 1 << endl;
			return;
		}

		CSv = CSn;

		for (int i = 0; i < ecoliSet[mj].nds.size(); i++) {
			ecoliSet[mi].nds.push_back(ecoliSet[mj].nds[i]);
		}

		vector<ecoli>::iterator imj = ecoliSet.begin();
		imj += mj;
		ecoliSet.erase(imj);
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
}

void ecoli::ecoli_cluster2() {
	vector<ecoli> ecoliSet;
	vector<vector<float>> dTable;
	ecoliSet = init_ecoli();
	for (int i = 0; i < ecoliSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < ecoliSet.size(); j++) {
			if (j > i)
				temp.push_back(ecoliDistance(ecoliSet[i], ecoliSet[j]));
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
		int mi = 0, mj = 0;
		vector<vector<float>> f;
		f.clear();

		for (int i = 0; i < ecoliSet.size(); i++) {
			vector<float> temp;
			for (int j = 0; j < ecoliSet.size(); j++) {
				if (j > i) {
					float sump = 0;
					for (vector<ecolinode>::iterator it1 = ecoliSet[i].nds.begin(); it1 < ecoliSet[i].nds.end(); it1++)
						for (vector<ecolinode>::iterator it2 = ecoliSet[j].nds.begin(); it2 < ecoliSet[j].nds.end(); it2++)
							sump += P[(*it1).id - 1][(*it2).id - 1];
					temp.push_back(sump / (ecoliSet[i].nds.size()*ecoliSet[j].nds.size()));
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
			cout << f.size() << endl;
			return;
		}

		if (maxf > 0.5&&ecoliSet.size() == 1) {
			cout << "error!" << endl;
			return;
		}

		for (vector<ecolinode>::iterator it1 = ecoliSet[mi].nds.begin(); it1 < ecoliSet[mi].nds.end(); it1++)
			for (vector<ecolinode>::iterator it2 = ecoliSet[mj].nds.begin(); it2 < ecoliSet[mj].nds.end(); it2++) {
				P[(*it1).id - 1][(*it2).id - 1] = 1 - P[(*it1).id - 1][(*it2).id - 1];
			}

		for (int i = 0; i < ecoliSet[mj].nds.size(); i++) {
			ecoliSet[mi].nds.push_back(ecoliSet[mj].nds[i]);
		}

		vector<ecoli>::iterator imj = ecoliSet.begin();
		imj += mj;
		ecoliSet.erase(imj);

	}
}

