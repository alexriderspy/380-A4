#include <bits/stdc++.h>
#include <omp.h>
#include "library.hpp"

using namespace std;

#define vi vector<int>
#define MAX_VAL INT_MAX

struct block{
    int row,col;
    vi matrix;
};

int main(int argc, char ** argv){
    unsigned char bytes[4];
    unsigned char bytes1[2];

    int n = 0, m = 0, k = 0, starti, startj;

    FILE *fp=fopen(argv[1],"rb");
    
    fread(bytes,4,1,fp);
    n  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    m  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    k  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    
    
    int cnt=0;
    map<pair<int,int>,vi>vv2;
    vector<block>vv(k);
    
    int i=0,j=0;
    while(cnt < k){
        fread(bytes, 4,1,fp);
        starti = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
        fread(bytes, 4,1,fp);
        startj = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));

        vi mat(m*m);

        for (int ti=0;ti<m;++ti){
            for(int tj=0;tj<m;++tj){
                fread(bytes1,2,1,fp);
                mat[ti*m+tj] = (bytes1[0] | (bytes1[1]<<8));       
            }
        }
        vv2[{starti, startj}] = mat;
        vv[cnt] = {starti, startj, mat};
        cnt++;
    }

    map<pair<int,int>, vi> mp;

    for(i=0;i<k;++i){
        for(j=0;j<k;++j){
            int r=-1, c=-1, index=-1;
            if(vv[i].row == vv[j].row && vv[i].col <= vv[j].col){
                r=vv[i].col, c=vv[j].col, index=vv[i].row; 
                if (mp.find({r,c})!=mp.end()){
                    mp[{r,c}].push_back(index);
                }else{
                    mp[{r,c}] = {index};
                }
            }
            if(vv[i].row == vv[j].col && vv[i].col <= vv[j].row){
                r=vv[i].col, c=vv[j].row, index=vv[i].row;
                if (mp.find({r,c})!=mp.end()){
                    mp[{r,c}].push_back(index);
                }else{
                    mp[{r,c}] = {index};
                }
            }
            if(vv[i].col == vv[j].row && vv[i].row <= vv[j].col){
                r=vv[i].row, c=vv[j].col, index = vv[i].col;
                if (mp.find({r,c})!=mp.end()){
                    mp[{r,c}].push_back(index);
                }else{
                    mp[{r,c}] = {index};
                }
            }
            if(vv[i].col == vv[j].col && vv[i].row <= vv[j].row){
                r=vv[i].row, c=vv[j].row, index=vv[i].col;
                if (mp.find({r,c})!=mp.end()){
                    mp[{r,c}].push_back(index);
                }else{
                    mp[{r,c}] = {index};
                }
            }
        }
    }

    int mpsize = mp.size();
    vector<block> finals(mpsize);
    int index=0;

    for(auto &x:mp){
        finals[index++] = {x.first.first,x.first.second,x.second};
    }

    vector<block> ans(mpsize);

    int total = 0;

    vector<int>is0(mpsize);

    #pragma omp parallel for num_threads(64)
    //{
    for(int i=0;i<mpsize;++i){
        #pragma omp task
        {
        vi out(m*m);
        for(int j=0;j<finals[i].matrix.size();++j){
            vi inn(m*m);
            vi fm,sm;

            if (finals[i].row <= finals[i].matrix[j] && finals[i].col <= finals[i].matrix[j]){
                fm = vv2[{finals[i].row,finals[i].matrix[j]}];
                sm = vv2[{finals[i].col,finals[i].matrix[j]}]; //fm.smT

                for(int ii=0;ii<m;ii++)
                    for(int kk=0;kk<m;kk++)
                        for(int jj=0;jj<m;jj++)
                            inn[ii*m + jj]=min(MAX_VAL,Outer(inn[ii*m+jj],min(MAX_VAL,Inner(fm[ii*m+kk],sm[kk+jj*m]))));

            }else if(finals[i].row <= finals[i].matrix[j] && finals[i].col > finals[i].matrix[j]){
                fm = vv2[{finals[i].row,finals[i].matrix[j]}];
                sm = vv2[{finals[i].matrix[j],finals[i].col}]; //fm.sm

                for(int ii=0;ii<m;ii++)
                    for(int kk=0;kk<m;kk++)
                        for(int jj=0;jj<m;jj++)
                            inn[ii*m + jj]=min(MAX_VAL,Outer(inn[ii*m+jj],min(MAX_VAL,Inner(fm[ii*m+kk],sm[kk*m+jj]))));

            }else if(finals[i].row > finals[i].matrix[j] && finals[i].col <= finals[i].matrix[j]){
                fm = vv2[{finals[i].matrix[j],finals[i].row}];
                sm = vv2[{finals[i].col,finals[i].matrix[j]}]; //fmTsmT

                for(int ii=0;ii<m;ii++)
                    for(int kk=0;kk<m;kk++)
                        for(int jj=0;jj<m;jj++)
                            inn[ii*m + jj]=min(MAX_VAL,Outer(inn[ii*m+jj],min(MAX_VAL,Inner(fm[ii+kk*m],sm[kk+jj*m]))));

            }else{
                fm = vv2[{finals[i].matrix[j],finals[i].row}];
                sm = vv2[{finals[i].matrix[j],finals[i].col}]; //fmTsm

                for(int ii=0;ii<m;ii++)
                    for(int kk=0;kk<m;kk++)
                        for(int jj=0;jj<m;jj++)
                            inn[ii*m + jj]=min(MAX_VAL,Outer(inn[ii*m+jj],min(MAX_VAL,Inner(fm[ii+kk*m],sm[kk*m+jj]))));

            }

            for (int ii=0;ii<m*m;++ii){
                out[ii] = min(MAX_VAL,Outer(out[ii],inn[ii]));
            }
        }
        ans[i] = {finals[i].row, finals[i].col,out};
        }
    }
    //}

    for(i=0;i<mpsize;++i){
        if (count(ans[i].matrix.begin(),ans[i].matrix.end(),0)!=m*m){
            total ++;
        }else{
            is0[i]=1;
        }
    }

    ofstream file(argv[2],ios::binary);
    file.write((char*)&n,4);
    file.write((char*)&m,4);
    file.write((char*)&total,4);
    for (i=0;i<mpsize;++i){
        if(is0[i]!=1){
            file.write((char*)&ans[i].row,4);
            file.write((char*)&ans[i].col,4);
            for (j=0;j<m*m;++j){
                file.write((char*)(&(ans[i].matrix[j])),4);
            }
        }
    }
    return 0;

}