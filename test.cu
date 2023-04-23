#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>

// function object that represents the loop body
struct loop_body
{
    __host__ __device__
    void operator()(int i) const
    {
        //dothis();
    }
};

int main()
{
    const int N = 1024;

    // create counting iterator for loop indices
    thrust::counting_iterator<int> indices(0);

    // use transform to apply loop body to each index in parallel
    thrust::for_each(indices, indices + N, loop_body());

    return 0;
}
