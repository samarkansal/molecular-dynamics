#include <iostream>
#include<bits/stdc++.h>
using namespace std;

int main(){
	int t;
	cin>>t;
	bool *ans=new bool[t];
	for (int i = 0; i < t; ++i)
	{
		ans[i]=false;
		int n;
		cin>>n;
		if(n<11){
			continue;
		}
		char str[100];
		cin>>str;
		for (int j = 0; j < n-10; ++j)
		{
			if(str[j]=='8'){
				ans[i]=true;
				break;
			}
		}
	}
	for (int i = 0; i < t; ++i)
	{
		if(ans[i]){
			cout<<"YES\n";
		}
		else{
			cout<<"NO\n";
		}
	}
	return 0;
}