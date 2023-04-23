#include <bits/stdc++.h>

using namespace std;

#define ll long long
#define vl vector<ll>
#define MAX_VAL ((1LL<<32)-1LL)

struct block{
    int row,col;
    vl matrix;
};

bool cmp(const block& b1, const block& b2){
    if (b1.row == b2.row){
        return b1.col < b2.col;
    }
    return b1.row < b2.row;
}

int main(int argc, char ** argv){
    unsigned char bytes[4];
    unsigned char bytes1[2];

    int n = 0, m = 0, k1 = 0, starti, startj;

    FILE *fp=fopen(argv[1],"rb");
    
    fread(bytes,4,1,fp);
    n  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    m  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    k1  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    
    
    int cnt=0;
    map<pair<int,int>,vl>vv2_1;
    vector<block>vv_1(k1);
    
    int i=0,j=0;
    while(cnt < k1){
        fread(bytes, 4,1,fp);
        starti = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
        fread(bytes, 4,1,fp);
        startj = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));

        vl mat(m*m);
        //cout<<starti<<' '<<startj<<'\n';
        for (int ti=0;ti<m;++ti){
            for(int tj=0;tj<m;++tj){
                fread(bytes1,2,1,fp);
                mat[ti*m+tj] = (bytes1[0] | (bytes1[1]<<8)); 
                cout<<mat[ti*m+tj]<<' ';      
            }
            cout<<'\n';
        }
        vv2_1[{starti, startj}] = mat;
        vv_1[cnt] = {starti, startj, mat};
        cnt++;
    }

    int k2 = 0;

    fp=fopen(argv[2],"rb");
    
    fread(bytes,4,1,fp);
    n  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    m  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    k2  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    
    
    cnt=0;
    map<pair<int,int>,vl>vv2_2;
    vector<block>vv_2(k2);
    
    i=0,j=0;
    while(cnt < k2){
        fread(bytes, 4,1,fp);
        starti = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
        fread(bytes, 4,1,fp);
        startj = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));

        vl mat(m*m);
        //cout<<starti<<' '<<startj<<'\n';
        for (int ti=0;ti<m;++ti){
            for(int tj=0;tj<m;++tj){
                fread(bytes1,2,1,fp);
                mat[ti*m+tj] = (bytes1[0] | (bytes1[1]<<8));   
                cout<<mat[ti*m+tj]<<' ';    
            }
            cout<<'\n';
        }
        vv2_2[{starti, startj}] = mat;
        vv_2[cnt] = {starti, startj, mat};
        cnt++;
    }
    //cout<<'\n';
    map<pair<int,int>, vl> mp;

    for(i=0;i<k1;++i){
        for(j=0;j<k2;++j){
            int r, c, index;
            if(vv_1[i].col == vv_2[j].row){
                r=vv_1[i].row, c=vv_2[j].col, index = vv_1[i].col;
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

    
    
    //{
    for(int i=0;i<mpsize;++i){
        vl out(m*m);
        
        for(int j=0;j<finals[i].matrix.size();++j){
            vl inn(m*m);
            vl fm,sm;

            fm = vv2_1[{finals[i].row,finals[i].matrix[j]}];
            sm = vv2_2[{finals[i].matrix[j],finals[i].col}];
            
            for(int ii=0;ii<m;ii++)
                for(int kk=0;kk<m;kk++)
                    for(int jj=0;jj<m;jj++)
                        inn[ii*m + jj]=min(MAX_VAL,inn[ii*m+jj]+min(MAX_VAL,fm[ii*m+kk]*sm[kk*m+jj]));

            for (int ii=0;ii<m*m;++ii){
                out[ii] = min(MAX_VAL,out[ii]+inn[ii]);
            }
        }
        ans[i] = {finals[i].row, finals[i].col,out};
    }
    //}

    for(i=0;i<mpsize;++i){
        if (count(ans[i].matrix.begin(),ans[i].matrix.end(),0)!=m*m){
            total ++;
        }else{
            is0[i]=1;
        }
    }

    sort(ans.begin(),ans.end(), cmp);
    ofstream file(argv[3],ios::binary);
    file.write((char*)&n,4);
    file.write((char*)&m,4);
    file.write((char*)&total,4);
    for (i=0;i<mpsize;++i){
        if(is0[i]!=1){
            file.write((char*)&ans[i].row,4);
            file.write((char*)&ans[i].col,4);
            //cout<<ans[i].row<<' '<<ans[i].col<<'\n';
            for (j=0;j<m*m;++j){
                file.write((char*)(&(ans[i].matrix[j])),4);
                //cout<<ans[i].matrix[j]<<' ';
            }
            //cout<<'\n';
        }
    }
    return 0;
}