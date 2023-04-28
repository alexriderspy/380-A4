#include <bits/stdc++.h>
using namespace std;

#define uint unsigned int
#define ll long long
#define vl vector<ll>
#define MAX_VAL ((1LL << 32) - 1LL)

__global__ void matrixMul(int *vA, int *vvA, int *vB, int *vvB, uint *vC, int n, int m)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < (n / m) && col < (n / m))
    {
        for (int k = 0; k < (n / m); ++k)
        {
            //blocks to be multiplied are (row,k) and (k,col)
            if(vA[row*(n/m) + k] != 0 && vB[k*(n/m) + col]!=0){
            for(int ii=0;ii<m;++ii){
                for(int kk=0;kk<m;++kk){
                    for(int jj=0;jj<m;++jj){
                        vC[(row*(n/m) + col)*m*m+ ii*m+jj] += vvA[vA[row*(n/m) + k] -1 + ii*m+kk]  * vvB[vB[k*(n/m) + col]-1 + kk*m+jj];
                    }
                }
            }
            }
        }
    }
}

int main(int argc, char **argv)
{
    unsigned char bytes[4];
    unsigned char bytes1[2];

    int n = 0, m = 0, k1 = 0, starti, startj;

    FILE *fp = fopen(argv[1], "rb");

    fread(bytes, 4, 1, fp);
    n = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
    fread(bytes, 4, 1, fp);
    m = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
    fread(bytes, 4, 1, fp);
    k1 = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));

    int *vA;
    int *vvA;

    cudaMallocManaged(&vA, (n / m) * (n / m) * sizeof(int));
    cudaMallocManaged(&vvA, k1 * m * m * sizeof(int));

    int cnt = 0;

    while (cnt < k1)
    {
        fread(bytes, 4, 1, fp);
        starti = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
        fread(bytes, 4, 1, fp);
        startj = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));

        for (int ti = 0; ti < m; ++ti)
        {
            for (int tj = 0; tj < m; ++tj)
            {
                fread(bytes1, 2, 1, fp);
                vvA[cnt * m * m + ti * m + tj] = (bytes1[0] | (bytes1[1] << 8));
            }
        }
        vA[starti * (n / m) + startj] = cnt * m*m+1;
        cnt++;
    }
    // cout << "m1\n";
    // for (int i = 0; i < (n/m); ++i)
    // {
    //     for (int j = 0; j < (n/m); ++j)
    //     {
    //         cout << vA[i * (n/m) + j] << ' ';
    //     }
    //     cout << '\n';
    // }
    // for (int i=0;i<k1*m*m;++i) cout<<vvA[i]<<' ';
    // cout<<'\n';
    fp = fopen(argv[2], "rb");

    fread(bytes, 4, 1, fp);
    n = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
    fread(bytes, 4, 1, fp);
    m = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
    fread(bytes, 4, 1, fp);
    int k2 = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));

    int *vB;
    int *vvB;

    cudaMallocManaged(&vB, (n / m) * (n / m) * sizeof(int));
    cudaMallocManaged(&vvB, k2 * m * m * sizeof(int));

    cnt = 0;

    while (cnt < k2)
    {
        fread(bytes, 4, 1, fp);
        starti = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
        fread(bytes, 4, 1, fp);
        startj = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));

        for (int ti = 0; ti < m; ++ti)
        {
            for (int tj = 0; tj < m; ++tj)
            {
                fread(bytes1, 2, 1, fp);
                vvB[cnt * m * m + ti * m + tj] = (bytes1[0] | (bytes1[1] << 8));
            }
        }
        vB[starti * (n / m) + startj] = cnt * m*m+1;
        cnt++;
    }

    uint *vC;

    //cudaError_t err1;

    cudaMallocManaged(&vC, n * n * sizeof(uint)); //(1<<30)

    //err1 = cudaPeekAtLastError();
    //cudaDeviceSynchronize();
    //printf("Got CUDA error ... %s \n", cudaGetErrorString(err1));

    int threads = 32;
    int blocks = ((n/m) + threads - 1) / threads;
    //cout<<blocks<<'\n';
    dim3 THREADS(threads, threads);
    dim3 BLOCKS(blocks, blocks);

    matrixMul<<<BLOCKS, THREADS>>>(vA,vvA,vB,vvB,vC,n,m);
    // cout << "m2\n";
    // for (int i = 0; i < n; ++i)
    // {
    //     for (int j = 0; j < n; ++j)
    //     {
    //         cout << b[i * n + j] << ' ';
    //     }
    //     cout << '\n';
    // }

    //err1 = cudaPeekAtLastError();
    cudaDeviceSynchronize();
    //printf("Got CUDA error ... %s \n", cudaGetErrorString(err1));


    int total = 0;
    vector<int> indices;

    for (int i  = 0; i< n*n; i+= m*m){
        int f=0;
        for(int j=i;j<m*m+i;++j){
            if(vC[j]!=0){
                f=1;break;
            }
        }
        if(f==1) {total++;indices.push_back(i);}
    }
    // cout << "m3\n";
    // for (int i = 0; i < n*n; i+=m*m)
    // {
    //     for (int j = i; j < m*m + i; ++j)
    //     {
    //         cout << vC[j] << ' ';
    //     }
    //     cout << '\n';
    // }

    ofstream file(argv[3], ios::binary);
    file.write((char *)&n, 4);
    file.write((char *)&m, 4);
    file.write((char *)&total, 4);
    //cout<<total<<'\n';
    for (int i = 0; i < indices.size(); ++i)
    {
        int i1 = indices[i]/(m*m);
        int r = i1/(n/m);
        int c=i1%(n/m);
        file.write((char *)&(r), 4);
        file.write((char *)&(c), 4);
        //cout<<r<<' '<<c<<'\n';
        // cout<<ans[i].row<<' '<<ans[i].col<<'\n';
        for (int k = i1 * m *m; k < i1*m*m+m*m; ++k)
        {
            file.write((char *)(&vC[k]), 4);
            
        }
    }
    return 0;

}