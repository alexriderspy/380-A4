#include <bits/stdc++.h>

using namespace std;

#define ll long long
#define vl vector<ll>
#define MAX_VAL ((1LL << 32) - 1LL)

__global__ void matrixMul(int *a, int *b, ll *c, int n)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    if (row < n && col < n)
    {
        ll tmp = 0;
        for (int i = 0; i < n; i++)
        {
            tmp = min(MAX_VAL,tmp +  (ll)a[row * n + i] * (ll)b[i * n + col]);
        }
        c[row * n + col] = tmp;
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

    int *a;

    cudaMallocManaged(&a, n * n); //(1<<30)
    int cnt = 0;

    while (cnt < k1)
    {
        fread(bytes, 4, 1, fp);
        starti = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
        fread(bytes, 4, 1, fp);
        startj = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
        // cout<<starti<<' '<<startj<<'\n';
        for (int ti = 0; ti < m; ++ti)
        {
            for (int tj = 0; tj < m; ++tj)
            {
                fread(bytes1, 2, 1, fp);
                // cout<<(bytes[0] | (bytes[1] << 8))<<" same as ";
                a[(starti * m + ti) * n + (startj * m + tj)] = (bytes1[0] | (bytes1[1] << 8));
                // cout<<a[((starti)*m+ti)*n+(startj*m+tj)]<<' ';
            }
            //   cout<<'\n';
        }
        // cout<<'\n';
        cnt++;
    }
    cout << "m1\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << a[i * n + j] << ' ';
        }
        cout << '\n';
    }
    fp = fopen(argv[2], "rb");

    fread(bytes, 4, 1, fp);
    n = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
    fread(bytes, 4, 1, fp);
    m = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
    fread(bytes, 4, 1, fp);
    int k2 = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));

    int *b;
    ll *c;

    cudaMallocManaged(&b, n * n); //(1<<30)
    cudaMallocManaged(&c, n * n); //(1<<30)
    cnt = 0;

    while (cnt < k2)
    {
        fread(bytes, 4, 1, fp);
        starti = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
        fread(bytes, 4, 1, fp);
        startj = (bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24));
        // cout<<starti<<' '<<startj<<'\n';
        for (int ti = 0; ti < m; ++ti)
        {
            for (int tj = 0; tj < m; ++tj)
            {
                fread(bytes1, 2, 1, fp);
                b[(starti * m + ti) * n + (startj * m + tj)] = (bytes1[0] | (bytes1[1] << 8));
                // cout<<b[((starti)*m+ti)*n+(startj*m+tj)]<<' ';
            }
            // cout<<'\n';
        }
        // cout<<'\n';
        cnt++;
    }
    cout << "m2\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << b[i * n + j] << ' ';
        }
        cout << '\n';
    }

    int threads = 32;
    int blocks = (n + threads - 1) / threads;

    dim3 THREADS(threads, threads);
    dim3 BLOCKS(blocks, blocks);

    matrixMul<<<BLOCKS, THREADS>>>(a, b, c, n);
    cudaDeviceSynchronize();

    int total = 0;
    vector<pair<int, int>> indices;

    for (int i = 0; i < (n / m); ++i)
    {
        for (int j = 0; j < (n / m); ++j)
        {
            int f = 0;
            for (int k = i * m; k < i * m + m; ++k)
            {
                for (int l = j * m; l < j * m + m; ++l)
                {
                    if (c[k * n + l] != 0)
                    {
                        f = 1;
                        break;
                    }
                }
                if (f == 1)
                    break;
            }
            if (f == 1)
            {
                total++;
                indices.push_back(make_pair(i, j));
                // cout<<i<<' '<<j<<'\n';
            }
        }
    }
    cout << "m3\n";
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            cout << c[i * n + j] << ' ';
        }
        cout << '\n';
    }

    ofstream file(argv[3], ios::binary);
    file.write((char *)&n, 4);
    file.write((char *)&m, 4);
    file.write((char *)&total, 4);
    for (int i = 0; i < indices.size(); ++i)
    {
        int i1 = indices[i].first;
        int i2 = indices[i].second;
        file.write((char *)&i1, 4);
        file.write((char *)&i2, 4);
        // cout<<ans[i].row<<' '<<ans[i].col<<'\n';
        for (int k = i1 * m; k < i1 * m + m; ++k)
        {
            for (int l = i2 * m; l < i2 * m + m; ++l)
            {
                file.write((char *)(&c[k * n + l]), 4);
            }
        }
    }
    return 0;
}
