#pragma once
#include<iostream>
using namespace std;
#include"mysql.h"

class dbconnect
{
public:
	dbconnect();
	~dbconnect();
	MYSQL *get_con() {
		return con;
	};
private:
	MYSQL * con;
};

dbconnect::dbconnect() {

	//database configuartion
	char dbuser[30] = "root";
	char dbpasswd[30] = "wp19970426";
	char dbip[30] = "localhost";
	char dbname[50] = "summer_test";

	con = mysql_init((MYSQL*)0);

	if (con != NULL && mysql_real_connect(con, dbip, dbuser, dbpasswd, dbname, 3306, NULL, 0)) {
		//cout << "连接成功！\n";
	}
	else {
		MessageBoxA(NULL, "Unable to connect the database,check your configuration!", "", NULL);
	}
}

dbconnect::~dbconnect() {
	mysql_close(con);
}

