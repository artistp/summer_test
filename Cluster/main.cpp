
#include"iris.h"
#include"wine.h"
#include"ecoli.h"
#include<iostream>
#include<vector>
#include<stdio.h>
using namespace std;

int main() {
	iris *ir = new iris();
	wine *wi = new wine();
	ecoli *ec = new ecoli();
	//ir->iris_cluster();
	//ir->iris_cluster2();

	//wi->wine_cluster();
	//wi->wine_cluster2();

	//ec->ecoli_cluster();
	ec->ecoli_cluster2();

	/*int i = 0;
	vector<iris>  ary;
	ary=ir->init_iris();
	for (i = 0; i < ary.size(); i++) {
		cout << ary[i].nds[0].id << "  ";
		cout << ary[i].nds[0].sepal_length << "  ";
		cout << ary[i].nds[0].sepal_width << "  ";
		cout << ary[i].nds[0].petal_length << "  ";
		cout << ary[i].nds[0].petal_width << "  ";
		cout << ary[i].nds[0].category << "  ";
		cout << endl;
	}*/
	system("PAUSE ");
	return 0;
}
