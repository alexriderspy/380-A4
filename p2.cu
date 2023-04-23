#include <stdio.h>

__global__ void add(int *a, int *b, int *c)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

int main()
{
    int a[3] = {1, 2, 3};
    int b[3] = {4, 5, 6};
    int c[3];

    int *dev_a, *dev_b, *dev_c;

    cudaMalloc((void**)&dev_a, 3 * sizeof(int));
    cudaMalloc((void**)&dev_b, 3 * sizeof(int));
    cudaMalloc((void**)&dev_c, 3 * sizeof(int));

    cudaMemcpy(dev_a, a, 3 * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, 3 * sizeof(int), cudaMemcpyHostToDevice);

    add<<<1,3>>>(dev_a, dev_b, dev_c);

    cudaMemcpy(c, dev_c, 3 * sizeof(int), cudaMemcpyDeviceToHost);

    for (int i = 0; i < 3; i++)
        printf("%d ", c[i]);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);

    return 0;
}
