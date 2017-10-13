/**
 * Sean McGlincy
 * Parallel Systems
 * Assignment 1
 *
 * Dependencies: I am using CLION as an IDE which uses CMAKE 3.8, and GCC, C+11
 * gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-11)
 * Running on Cento 7
 *
 *
 * Program:

 *
 * */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <assert.h>

using namespace std;
using namespace chrono;

#ifdef _OPENMP
#include <omp.h>
#endif

void get_rank_thread_count(int *my_rank, int *thread_count) {
#ifdef _OPENMP
    *my_rank = omp_get_thread_num();
    *thread_count = omp_get_num_threads();
#else
    *my_rank = 1;
    *thread_count = 0;
#endif
}


void clock(high_resolution_clock::time_point *array, int* time_samples){
    for (int i = 0; i < *time_samples; i++) {
        array[i] = high_resolution_clock::now();
    }
}

double calculate_time(high_resolution_clock::time_point *start, high_resolution_clock::time_point *end, int *time_samples){
    // Average time and convert to Micro Sec; 1 sec = 1,000,000 micro sec
    double total = 0;
    for (int i = 0; i < *time_samples; i++) {
        chrono::duration<double, std::milli> diff = end[i] - start[i];  // Time in Micro Sec
        total += diff.count();
    }
    return total / *time_samples;
}

void print(int sample_size, int min, int max, int bucketSize, vector<int> *data) {
    //////// Print //////////////
    cout << "Sample Size: " << sample_size << endl;
    cout << "Min Value: " << min << endl;
    cout << "Max Value: " << max << endl;
    cout << "bucket Range: " << bucketSize - 1 << endl;
    for (int i = 0; i < data->size(); ++i) {
        cout << "[" << min + (i * bucketSize) << ", " << min + ((i + 1) * bucketSize) - 1 << "] : " << data->at(i)
             << endl;
    }
}



int check_user_number(char *argv){
    char *endptr;
    int intervalSize = strtol(argv, &endptr, 10);
    if (!*argv || *endptr)
        cerr << "Invalid number " << argv << '\n';
    return intervalSize;
}


int main(int argc, char *argv[]) {
    cout << "Starting Program" << endl;


    //////// Start Clock //////////////
    // Use Chrono for high grade clock
    int time_samples = 5;
    high_resolution_clock::time_point clock_start[time_samples];
    high_resolution_clock::time_point clock_end[time_samples];
    clock(clock_start, &time_samples);


    //////// USER INPUT//////////////
    if (argc != 4) {
        cout << "Error Error" << endl;
        cout << "Please provide: binary data file and interval size" << endl;
        exit(1);
    }

    assert(argc == 4);
    string filePath = argv[1]; //Filename

    int intervalSize =  check_user_number(argv[2]);  // Interval Size
    assert(intervalSize > 0 );

    int numThreads =  check_user_number(argv[3]);  // NumThreads
    assert(numThreads > 0 );


    //////// Variables //////////////
    int unit = sizeof(int);
    int bufferSize = unit * 10000;  // read in chuncks of file
    int size = bufferSize;
    int bucketSize = 0;
    int fileLength = 0;
    int max = numeric_limits<int>::min();
    int min = numeric_limits<int>::max();
    vector<int> buffer(bufferSize, 0);
    vector<int> data(intervalSize, 0);
    int my_rank = 0;
    int thread_count = 1;



    //////// OPEN FILE //////////////
    ifstream fileInput;
    fileInput.open(filePath, ios::binary);
    if (fileInput.is_open()) {

        //////// Read FILE //////////////
        // Get file size
        fileInput.seekg(0, ios::end);
        fileLength = fileInput.tellg();
        fileInput.seekg(0);

        //////// First Pass:  Min Max   //////////////
        for (int i = 0; i < fileLength; i += bufferSize) {
            //  Check if buffer is less then remainder of file
            if (fileLength - i < bufferSize) {
                size = fileLength - i;
            }

            fileInput.read((char *) &buffer[0], size);
# pragma omp parallel num_threads(numThreads) default(none) private(i, my_rank, thread_count) shared(buffer, min, max, size, unit)
            {
                get_rank_thread_count(&my_rank, &thread_count);
                int local_min = numeric_limits<int>::max();
                int local_max = numeric_limits<int>::min();
// Find max and min
# pragma omp for schedule(static, ( size / unit ) / thread_count )
                for (int i = 0; i < size / unit; ++i) {
                    if (local_min > buffer[i])
                        local_min = buffer[i];
                    if (local_max < buffer[i])
                        local_max = buffer[i];
                }
# pragma omp critical
                {
                    if (min > local_min)
                        min = local_min;
                    if (max < local_max)
                        max = local_max;
                }
            }

        }


//    //  Start OMP
//    #pragma omp parallel num_threads(numThreads) default(none) private(i) shared(buffer, size, unit, min, max) reduction(max: max) reduction(min: min){
//
//
//        //////// First Pass:  Min Max   //////////////
//        for(int i = 0; i < fileLength; i += bufferSize){
//            //  Check if buffer is less then remainder of file
//            if(fileLength - i < bufferSize ) {
//               size = fileLength - i;
//            }
//
//            fileInput.read((char *) &buffer[0], size);
//
//            }
//
//
//
//        }



        // Reset
        fileInput.seekg(0);
        size = bufferSize;

        // Check Values and calculte bucket size
        int range = abs(max - min);
        if (intervalSize > range){  // if user requesting to many buckets
            intervalSize = range;
            data.resize(intervalSize);
        }
        bucketSize =  range / intervalSize;
        bucketSize++;

        //////// Second Pass:  add the values   //////////////
        for(int i = 0; i < fileLength; i += bufferSize){
            if(fileLength - i < bufferSize ) {   //  Check if buffer is less then remainder of file
                size = fileLength - i;
            }

            fileInput.read((char *) &buffer[0], size);

            /**
             * Second attempt at non blocking critical section
             * Algorithm from: https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
             * Algorithm seems slower then using a critical section and is more dangourous and harder to read
             */

/*
            int* local_data;
# pragma omp parallel num_threads(numThreads) default(none) private(my_rank)  shared(data, local_data, buffer, size, unit, intervalSize, min, bucketSize, thread_count)
            {
                get_rank_thread_count(&my_rank, &thread_count);

# pragma omp single
                {
                    local_data = new int[thread_count * intervalSize];
                    for(int j = 0; j< thread_count * intervalSize; j++)
                        local_data[j] = 0;


                }
                // Process Data
# pragma omp  for schedule(static, ( size / unit ) / thread_count )

                for (int i = 0; i < size / unit; i++) {
                    local_data[(my_rank * intervalSize)  + ((buffer[i] - min) / bucketSize)]++;
                }
# pragma omp for schedule(static, intervalSize / thread_count )
                for(int i = 0; i< intervalSize; i++){
                    for (int process = 0; process < thread_count; process++) {
                        data[i] += local_data[intervalSize * process + i];
                    }
                }

            } // End parrallel
            delete[] local_data;
        }
        */


/* Original Reduction */
# pragma omp parallel num_threads(numThreads) default(none) private(my_rank, thread_count)  shared(data, buffer, size, unit, intervalSize, min, bucketSize)
            {
                // Process Data
                vector<int> local_data(intervalSize, 0);
                get_rank_thread_count(&my_rank, &thread_count);

# pragma omp  for schedule(static, ( size / unit ) / thread_count )
                for (int i = 0; i < size / unit; i++) {
                    local_data[(buffer[i] - min) / bucketSize]++;
                }
# pragma omp critical
                {
                    for (int i = 0; i < intervalSize; i++) {
                        data[i] += local_data[i];
                    }
                }

            }
        }

        fileInput.close();
        //////// End FILE //////////////

    } else {
        cout << "Can Not open file..." << endl;
        exit(1);
    }



    //////// Print //////////////
    print(fileLength / unit, min, max, bucketSize, &data);

    // Release resources
    buffer.clear();
    data.clear();


    ////////  END CLOCK //////////////
    //////// GET TIME //////////////
    clock(clock_end, &time_samples);
    double total_time =  calculate_time(clock_start, clock_end, &time_samples);
    cout << "AVG Time: " << total_time  << " Milli Seconds" << endl;
    cout << "Program complete!" << endl;

    return 0;
}


