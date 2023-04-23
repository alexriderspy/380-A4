#include <bits/stdc++.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>

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

__global__ void matrixMul(){

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
    
    
    int blks = (n/m);

    int *v1_ind;
    int *v1_val;
    cudaMallocManaged(&v1_ind, blks*blks); //(1<<30)
    cudaMallocManaged(&v1_val, k1*8*(m*m)); // row, col, m*m 

    int step = 8*m*m;
    int cnt=0;
    // map<pair<int,int>,vl>vv2_1;
    // vector<block>vv_1(k1);
    
    int i=0,j=0;
    while(cnt < k1){
        fread(bytes, 4,1,fp);
        starti = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
        fread(bytes, 4,1,fp);
        startj = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));

        v1_ind[starti*blks+startj] = cnt;
        int s = cnt;
        v1_val[s++] = starti;
        v1_val[s++] = startj;

        for (int ti=0;ti<m;++ti){
            for(int tj=0;tj<m;++tj){
                fread(bytes1,2,1,fp);
                v1_val[s++] = (bytes1[0] | (bytes1[1]<<8)); 
            }
        }
        //vv2_1[{starti, startj}] = mat;
        //vv_1[cnt] = {starti, startj, mat};
        cnt++;
    }

    fp=fopen(argv[2],"rb");
    
    fread(bytes,4,1,fp);
    n  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    m  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
    fread(bytes,4,1,fp);
    int k2  = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));

    int *v2_ind;
    int *v2_val;
    cudaMallocManaged(&v2_ind, blks*blks); //(1<<30)
    cudaMallocManaged(&v2_val, k2*8*(m*m)); // row, col, m*m 

    step = 8*m*m;
    cnt=0;
    // map<pair<int,int>,vl>vv2_1;
    // vector<block>vv_1(k1);
    
    i=0,j=0;
    while(cnt < k1){
        fread(bytes, 4,1,fp);
        starti = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));
        fread(bytes, 4,1,fp);
        startj = (bytes[0] | (bytes[1]<<8) | (bytes[2]<<16) | (bytes[3]<<24));

        v2_ind[starti*blks+startj] = cnt;
        int s = cnt;
        v2_val[s++] = starti;
        v2_val[s++] = startj;

        for (int ti=0;ti<m;++ti){
            for(int tj=0;tj<m;++tj){
                fread(bytes1,2,1,fp);
                v2_val[s++] = (bytes1[0] | (bytes1[1]<<8)); 
            }
        }
        //vv2_1[{starti, startj}] = mat;
        //vv_1[cnt] = {starti, startj, mat};
        cnt++;
    }

    // cout<<'\n';
    map<pair<int,int>, vl> mp;

    int total_mp = 0;
    //do -1 before using 
    for(i=0;i<k1*step;i+=step){
        for(j=0;j<k2*step;j+=step){
            int r, c, index;
            if(v1_val[i+1] == v2_val[j]){
                r=v1_val[i], c=v2_val[j+1], index = v1_val[i+1];

                if (mp.find({r,c})!=mp.end()){
                    mp[{r,c}].push_back(index);
                    total_mp+=3;
                }else{
                    mp[{r,c}] = {index};
                    total_mp++;
                }
            }
        }
    }

    int mpsize = mp.size();

    int *f_ind;
    int *f_val;

    cudaMallocManaged(&f_ind, mpsize+1); //+1 to have the last entry
    cudaMallocManaged(&f_val, total_mp);

    //vector<block> finals(mpsize);
    int index_val=0, index_ind =0;

    for(auto &x:mp){
        f_ind[index_ind++] = index_val;
        f_val[index_val++] = x.first.first;
        f_val[index_val++] = x.first.second;
        for (int i=0;i<x.second.size();++i)
            f_val[index_val++] = x.second[i];
    }

    f_ind[index_ind] = index_val;

    vector<block> ans(mpsize);

    int total = 0;

    vector<int>is0(mpsize);

    int THREADS = 1000;
    int BLOCKS = (mpsize + THREADS -1)/THREADS;


    matrixMul<<<BLOCKS, THREADS>>>();
    
    __global__ void matrixMul();{

        int tid = blockIdx.x * blockDim.x + threadIdx.x;

        if (tid < mpsize){
            vl out(m*m);
            
            int r=f_val[f_ind[tid]], c=f_val[f_ind[tid]+1];
            for(int j=f_ind[tid]+2;j<f_ind[tid+1];++j){

                int index = f_val[j];
                vl inn(m*m);
                vl fm,sm;
                

                //multiply (r,index) * (index, c)
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
        
    }
    // //{
    // for(int i=0;i<mpsize;++i){
    //     vl out(m*m);
        
    //     for(int j=0;j<finals[i].matrix.size();++j){
    //         vl inn(m*m);
    //         vl fm,sm;

    //         fm = vv2_1[{finals[i].row,finals[i].matrix[j]}];
    //         sm = vv2_2[{finals[i].matrix[j],finals[i].col}];
            
    //         for(int ii=0;ii<m;ii++)
    //             for(int kk=0;kk<m;kk++)
    //                 for(int jj=0;jj<m;jj++)
    //                     inn[ii*m + jj]=min(MAX_VAL,inn[ii*m+jj]+min(MAX_VAL,fm[ii*m+kk]*sm[kk*m+jj]));

    //         for (int ii=0;ii<m*m;++ii){
    //             out[ii] = min(MAX_VAL,out[ii]+inn[ii]);
    //         }
    //     }
    //     ans[i] = {finals[i].row, finals[i].col,out};
    // }
    // //}

    // for(i=0;i<mpsize;++i){
    //     if (count(ans[i].matrix.begin(),ans[i].matrix.end(),0)!=m*m){
    //         total ++;
    //     }else{
    //         is0[i]=1;
    //     }
    // }
    
    // thrust::sort(thrust::host, ans.begin(), ans.end(), cmp);
    // //sort(ans.begin(),ans.end(), cmp);
    // ofstream file(argv[3],ios::binary);
    // file.write((char*)&n,4);
    // file.write((char*)&m,4);
    // file.write((char*)&total,4);
    // for (i=0;i<mpsize;++i){
    //     if(is0[i]!=1){
    //         file.write((char*)&ans[i].row,4);
    //         file.write((char*)&ans[i].col,4);
    //         //cout<<ans[i].row<<' '<<ans[i].col<<'\n';
    //         for (j=0;j<m*m;++j){
    //             file.write((char*)(&(ans[i].matrix[j])),4);
    //             //cout<<ans[i].matrix[j]<<' ';
    //         }
    //         //cout<<'\n';
    //     }
    // }
    return 0;
}