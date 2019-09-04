#pragma once
#include"mysql.h"
#include<iostream>
#include<vector>
#include<stack>
#include<cmath>
#include<string>
using namespace std;
#include"dbconnect.h"
class winenode {
public:
	int id;
	int category;
	float Alcohol;
	float Malic_acid;
	float Ash;
	float Alcalinity_of_ash;
	float Magnesium;
	float Total_phenols;
	float Flavanoids;
	float Nonflavanoid_phenols;
	float Proanthocyanins;
	float Color_intensity;
	float Hue;
	float diluted_wines;
	float Proline;
	winenode(int i,int categor,float Alcoho,float Malic_aci,float As,float Alcalinity_of_as,float Magnesiu,float Total_phenol,
	float Flavanoid,float Nonflavanoid_phenol,float Proanthocyanin,float Color_intensit,float Hu,float diluted_wine,float Prolin) {
		id=i;
		category=categor;
		Alcohol= Alcoho;
		Malic_acid= Malic_aci;
		Ash= As;
		Alcalinity_of_ash=Alcalinity_of_as;
		Magnesium= Magnesiu;
		Total_phenols=Total_phenol;
		Flavanoids=Flavanoid;
		Nonflavanoid_phenols=Nonflavanoid_phenol;
		Proanthocyanins=Proanthocyanin;
		Color_intensity=Color_intensit;
		Hue=Hu;
		diluted_wines=diluted_wine;
		Proline=Prolin;
	}
};
float AlcoholMax = INT_MIN;
float AlcoholMin = INT_MAX;
float Malic_acidMax = INT_MIN;
float Malic_acidMin = INT_MAX;
float AshMax = INT_MIN;
float AshMin = INT_MAX;
float Alcalinity_of_ashMax = INT_MIN;
float Alcalinity_of_ashMin = INT_MAX;
float MagnesiumMax = INT_MIN;
float MagnesiumMin = INT_MAX;
float Total_phenolsMax = INT_MIN;
float Total_phenolsMin = INT_MAX;
float FlavanoidsMax = INT_MIN;
float FlavanoidsMin = INT_MAX;
float Nonflavanoid_phenolsMax = INT_MIN;
float Nonflavanoid_phenolsMin = INT_MAX;
float ProanthocyaninsMax = INT_MIN;
float ProanthocyaninsMin = INT_MAX;
float Color_intensityMax = INT_MIN;
float Color_intensityMin = INT_MAX;
float HueMax = INT_MIN;
float HueMin = INT_MAX;
float diluted_winesMax = INT_MIN;
float diluted_winesMin = INT_MAX;
float ProlineMax = INT_MIN;
float ProlineMin = INT_MAX;
class wine
{
public:
	int id;
	int category;
	float Alcohol;
	float Malic_acid;
	float Ash;
	float Alcalinity_of_ash;
	float Magnesium;
	float Total_phenols;
	float Flavanoids;
	float Nonflavanoid_phenols;
	float Proanthocyanins;
	float Color_intensity;
	float Hue;
	float diluted_wines;
	float Proline;
	vector<winenode> nds;
	wine(int i,int categor,float Alcoho,float Malic_aci,float As,float Alcalinity_of_as,float Magnesiu,float Total_phenol,
		float Flavanoid,float Nonflavanoid_phenol,float Proanthocyanin,float Color_intensit,float Hu,float diluted_wine,float Prolin);
	wine();
	~wine();
	vector<wine> init_wine();
	float wineDistance(wine w1, wine w2);
	void wine_cluster();
	void wine_cluster2();
};

wine::wine() {

}
wine::~wine() {

}
wine::wine(int i,int categor,float Alcoho,float Malic_aci,float As,float Alcalinity_of_as,float Magnesiu,float Total_phenol,
	float Flavanoid,float Nonflavanoid_phenol,float Proanthocyanin,float Color_intensit,float Hu,float diluted_wine,float Prolin) {
	id = i;
	category = categor;
	Alcohol = Alcoho;
	Malic_acid = Malic_aci;
	Ash = As;
	Alcalinity_of_ash = Alcalinity_of_as;
	Magnesium = Magnesiu;
	Total_phenols = Total_phenol;
	Flavanoids = Flavanoid;
	Nonflavanoid_phenols = Nonflavanoid_phenol;
	Proanthocyanins = Proanthocyanin;
	Color_intensity = Color_intensit;
	Hue = Hu;
	diluted_wines = diluted_wine;
	Proline = Prolin;
	winenode t(id,category ,Alcohol ,Malic_acid ,Ash ,Alcalinity_of_ash ,Magnesium ,Total_phenols ,Flavanoids ,Nonflavanoid_phenols ,Proanthocyanins ,Color_intensity ,Hue ,diluted_wines ,Proline );
	nds.push_back(t);
}

float wine::wineDistance(wine w1, wine w2) {
	return (abs(w1.Alcohol - w2.Alcohol) / (AlcoholMax - AlcoholMin) +
		abs(w1.Malic_acid - w2.Malic_acid) / (Malic_acidMax - Malic_acidMin) +
		abs(w1.Ash - w2.Ash) / (AshMax - AshMin) +
		abs(w1.Alcalinity_of_ash - w2.Alcalinity_of_ash) / (Alcalinity_of_ashMax - Alcalinity_of_ashMin)+
		abs(w1.Magnesium - w2.Magnesium) / (MagnesiumMax - MagnesiumMin)+
		abs(w1.Total_phenols - w2.Total_phenols) / (Total_phenolsMax - Total_phenolsMin)+
		abs(w1.Flavanoids - w2.Flavanoids) / (FlavanoidsMax - FlavanoidsMin)+
		abs(w1.Nonflavanoid_phenols - w2.Nonflavanoid_phenols) / (Nonflavanoid_phenolsMax - Nonflavanoid_phenolsMin)+
		abs(w1.Proanthocyanins - w2.Proanthocyanins) / (ProanthocyaninsMax - ProanthocyaninsMin)+
		abs(w1.Color_intensity - w2.Color_intensity) / (Color_intensityMax - Color_intensityMin)+
		abs(w1.Hue - w2.Hue) / (HueMax - HueMin)+
		abs(w1.diluted_wines - w2.diluted_wines) / (diluted_winesMax - diluted_winesMin)+
		abs(w1.Proline - w2.Proline) / (ProlineMax - ProlineMin)) / 13;
}

vector<wine> wine::init_wine() {
	dbconnect *dbcon = new dbconnect();
	MYSQL * con = dbcon->get_con();
	MYSQL_RES *res;
	MYSQL_ROW row;
	char tmp[400];
	char tablename[50] = "wine";
	char * query = NULL;
	int rt;//return value;

	vector<wine> ary;

	sprintf_s(tmp, "select * from %s ", tablename);
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
		wine crount(atof(row[0]),atof(row[1]), atof(row[2]), atof(row[3]), atof(row[4]), atof(row[5]), atof(row[6]), atof(row[7]), atof(row[8]), atof(row[9]), atof(row[10]), atof(row[11]), atof(row[12]), atof(row[13]), atof(row[14]));
		if (AlcoholMax < atof(row[2]))
			AlcoholMax = atof(row[2]);
		if (AlcoholMin > atof(row[2]))
			AlcoholMin = atof(row[2]);
		if (Malic_acidMax < atof(row[3]))
			Malic_acidMax = atof(row[3]);
		if (Malic_acidMin > atof(row[3]))
			Malic_acidMin = atof(row[3]);
		if (AshMax < atof(row[4]))
			AshMax = atof(row[4]);
		if (AshMin > atof(row[4]))
			AshMin = atof(row[4]);
		if (Alcalinity_of_ashMax < atof(row[5]))
			Alcalinity_of_ashMax = atof(row[5]);
		if (Alcalinity_of_ashMin > atof(row[5]))
			Alcalinity_of_ashMin = atof(row[5]);
		if (MagnesiumMax < atof(row[6]))
			MagnesiumMax = atof(row[6]);
		if (MagnesiumMin > atof(row[6]))
			MagnesiumMin = atof(row[6]);
		if (Total_phenolsMax < atof(row[7]))
			Total_phenolsMax = atof(row[7]);
		if (Total_phenolsMin > atof(row[7]))
			Total_phenolsMin = atof(row[7]);
		if (FlavanoidsMax < atof(row[8]))
			FlavanoidsMax = atof(row[8]);
		if (FlavanoidsMin > atof(row[8]))
			FlavanoidsMin = atof(row[8]);
		if (Nonflavanoid_phenolsMax < atof(row[9]))
			Nonflavanoid_phenolsMax = atof(row[9]);
		if (Nonflavanoid_phenolsMin > atof(row[9]))
			Nonflavanoid_phenolsMin = atof(row[9]);
		if (ProanthocyaninsMax < atof(row[10]))
			ProanthocyaninsMax = atof(row[10]);
		if (ProanthocyaninsMin > atof(row[10]))
			ProanthocyaninsMin = atof(row[10]);
		if (Color_intensityMax < atof(row[11]))
			Color_intensityMax = atof(row[11]);
		if (Color_intensityMin > atof(row[11]))
			Color_intensityMin = atof(row[11]);
		if (HueMax < atof(row[12]))
			HueMax = atof(row[12]);
		if (HueMin > atof(row[12]))
			HueMin = atof(row[12]);
		if (diluted_winesMax < atof(row[13]))
			diluted_winesMax = atof(row[13]);
		if (diluted_winesMin > atof(row[13]))
			diluted_winesMin = atof(row[13]);
		if (ProlineMax < atof(row[14]))
			ProlineMax = atof(row[14]);
		if (ProlineMin > atof(row[14]))
			ProlineMin = atof(row[14]);
		ary.push_back(crount);
	}
	mysql_free_result(res);
	free(dbcon);
	return ary;
}

void wine::wine_cluster() {
	vector<wine> wineSet;
	vector<vector<float>> dTable;
	wineSet = init_wine();
	for (int i = 0; i < wineSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < wineSet.size(); j++) {
			if (j > i)
				temp.push_back(wineDistance(wineSet[i], wineSet[j]));
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
			cout << wineSet.size() << endl;
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

		if (CSv > CSn&&wineSet.size() == 1) {
			cout << 1 << endl;
			return;
		}

		if (CSv <= CSn) {
			cout << wineSet.size() + 1 << endl;
			return;
		}

		CSv = CSn;

		for (int i = 0; i < wineSet[mj].nds.size(); i++) {
			wineSet[mi].nds.push_back(wineSet[mj].nds[i]);
		}

		vector<wine>::iterator imj = wineSet.begin();
		imj += mj;
		wineSet.erase(imj);
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

void wine::wine_cluster2() {
	vector<wine> wineSet;
	vector<vector<float>> dTable;
	wineSet = init_wine();
	for (int i = 0; i < wineSet.size(); i++) {
		vector<float> temp;
		for (int j = 0; j < wineSet.size(); j++) {
			if (j > i)
				temp.push_back(wineDistance(wineSet[i], wineSet[j]));
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

		for (int i = 0; i < wineSet.size(); i++) {
			vector<float> temp;
			for (int j = 0; j < wineSet.size(); j++) {
				if (j > i) {
					float sump = 0;
					for (vector<winenode>::iterator it1 = wineSet[i].nds.begin(); it1 < wineSet[i].nds.end(); it1++)
						for (vector<winenode>::iterator it2 = wineSet[j].nds.begin(); it2 < wineSet[j].nds.end(); it2++)
							sump += P[(*it1).id - 1][(*it2).id - 1];
					temp.push_back(sump / (wineSet[i].nds.size()*wineSet[j].nds.size()));
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

		if (maxf > 0.5&&wineSet.size() == 1) {
			cout << "error!" << endl;
			return;
		}

		for (vector<winenode>::iterator it1 = wineSet[mi].nds.begin(); it1 < wineSet[mi].nds.end(); it1++)
			for (vector<winenode>::iterator it2 = wineSet[mj].nds.begin(); it2 < wineSet[mj].nds.end(); it2++) {
				P[(*it1).id - 1][(*it2).id - 1] = 1 - P[(*it1).id - 1][(*it2).id - 1];
			}

		for (int i = 0; i < wineSet[mj].nds.size(); i++) {
			wineSet[mi].nds.push_back(wineSet[mj].nds[i]);
		}

		vector<wine>::iterator imj = wineSet.begin();
		imj += mj;
		wineSet.erase(imj);
	}
}


