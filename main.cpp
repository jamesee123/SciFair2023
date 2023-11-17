#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <thread>

#define arraySize 100000

using std::string;

void b(int l) {
    std::cout << l << "\n";
}

int partition(int* arr, int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            std::swap(arr[i], arr[j]);
        }
    }

    std::swap(arr[i + 1], arr[high]);
    return (i + 1);
}

inline void selection(int* array)
{
    for (int i = 0; i < arraySize - 1; i++) {
        // Find the minimum element in the unsorted part of the array
        int minIndex = i;
        for (int j = i + 1; j < arraySize; j++) {
            if (array[j] < array[minIndex]) {
                minIndex = j;
            }
        }

        // Swap the found minimum element with the element at position i
        if (minIndex != i) {
            std::swap(array[i], array[minIndex]);
        }
    }
}

void quickR(int* arr, int low, int high)
{
    if (low < high) {
        int pivot = partition(arr, low, high);
        quickR(arr, low, pivot - 1);
        quickR(arr, pivot + 1, high);
    }
}

inline void quick(int* arr)
{
    int pivot = partition(arr, 0, 99);
    quickR(arr, 0, pivot - 1);
    quickR(arr, pivot + 1, 99);
}

void heapify(int* arr, int n, int i) {
    int largest = i; // Initialize the largest as the root
    int left = 2 * i + 1; // Left child
    int right = 2 * i + 2; // Right child

    // If the left child is larger than the root
    if (left < n && arr[left] > arr[largest])
        largest = left;

    // If the right child is larger than the largest so far
    if (right < n && arr[right] > arr[largest])
        largest = right;

    // If the largest is not the root
    if (largest != i) {
        std::swap(arr[i], arr[largest]);

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

inline void heap(int* arr)
{
    // Build a max heap
    for (int i = arraySize / 2 - 1; i >= 0; i--)
        heapify(arr, arraySize, i);

    // Extract elements one by one from the heap
    for (int i = arraySize - 1; i > 0; i--) {
        std::swap(arr[0], arr[i]); // Move the current root to the end
        heapify(arr, i, 0); // Call max heapify on the reduced heap
    }
}

inline void insertion(int* array)
{
    for (int i = 1; i < arraySize; i++) {
        int key = array[i];
        int j = i - 1;

        while (j >= 0 && array[j] > key) {
            array[j + 1] = array[j];
            j--;
        }

        array[j + 1] = key;
    }
}

inline void bubble(int* array)
{
    for (int i = 0; i < arraySize - 1; i++) {
        for (int j = 0; j < arraySize - i - 1; j++) {
            if (array[j] > array[j + 1]) {
                std::swap(array[j], array[j + 1]);
            }
        }
    }
}

void merge(int* arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    // Create temporary arrays
    std::vector<int> L(n1);
    std::vector<int> R(n2);

    // Copy data to temporary arrays L[] and R[]
    for (int i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (int j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    // Merge the temporary arrays back into arr[l..r]
    int i = 0;
    int j = 0;
    int k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    // Copy the remaining elements of L[], if there are any
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    // Copy the remaining elements of R[], if there are any
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

// Merge Sort function with a single parameter

// Main merge sort function
void mergeSortHelper(int* arr, int l, int r) {
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for large l and r
        int m = l + (r - l) / 2;

        // Sort first and second halves
        mergeSortHelper(arr, l, m);
        mergeSortHelper(arr, m + 1, r);

        // Merge the sorted halves
        merge(arr, l, m, r);
    }
}
inline void mergeSort(int* arr) {
    mergeSortHelper(arr, 0, arraySize - 1);
}

inline void counting(int* arr) {
    int max = arr[0];
    for (int i = 1; i < arraySize; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }

    int min = arr[0];
    for (int i = 1; i < arraySize; i++) {
        if (arr[i] < min) {
            min = arr[i];
        }
    }

    int range = max - min + 1;
    int* count = new int[range];
    int* output = new int[arraySize];

    for (int i = 0; i < range; i++) {
        count[i] = 0;
    }

    for (int i = 0; i < arraySize; i++) {
        count[arr[i] - min]++;
    }

    for (int i = 1; i < range; i++) {
        count[i] += count[i - 1];
    }

    for (int i = arraySize - 1; i >= 0; i--) {
        output[count[arr[i] - min] - 1] = arr[i];
        count[arr[i] - min]--;
    }

    for (int i = 0; i < arraySize; i++) {
        arr[i] = output[i];
    }

    delete[] count;
    delete[] output;
}

struct Bucket {
    std::vector<int> values;
};

inline void bucket(int* arr) {
    // Find the maximum and minimum values in the array
    int maxVal = arr[0];
    int minVal = arr[0];
    for (int i = 1; i < arraySize; i++) {
        if (arr[i] > maxVal) {
            maxVal = arr[i];
        }
        if (arr[i] < minVal) {
            minVal = arr[i];
        }
    }

    // Create an array of buckets
    int numBuckets = maxVal - minVal + 1;
    std::vector<Bucket> buckets(numBuckets);

    // Distribute the elements into buckets
    for (int i = 0; i < arraySize; i++) {
        int index = arr[i] - minVal;
        buckets[index].values.push_back(arr[i]);
    }

    // Sort each bucket and update the original array
    int k = 0;
    for (int i = 0; i < numBuckets; i++) {
        std::sort(buckets[i].values.begin(), buckets[i].values.end());
        for (int j = 0; j < buckets[i].values.size(); j++) {
            arr[k] = buckets[i].values[j];
            k++;
        }
    }
}


void runAlgorithm(int* arr, string algorithm, bool outCSV)
{
    std::chrono::milliseconds delay(1);
    if (!outCSV) {
        std::cout << algorithm << " Sort: " << std::flush;
    }
    auto start = std::chrono::high_resolution_clock::now();
    switch (algorithm[1]+algorithm[2]-algorithm[3])
    {
        case 108:
            selection(arr);
            break;
        case 123:
            quick(arr);
            break;
        case 86:
            heap(arr);
            break;
        case 124:
            insertion(arr);
            break;
        case 117:
            bubble(arr);
            break;
        case 112:
            mergeSort(arr);
            break;
        case 118:
            counting(arr);
            break;
        case 109:
            bucket(arr);
            break;
    }
    std::this_thread::sleep_for(delay);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::microseconds duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (!outCSV) {
        std::cout << duration.count() << std::endl;
    } else {
        std::cout << "," << duration.count();
    }
}

/*
Selection: 108
Quick:     123
Heap:      86
Insertion: 124
Bubble:    117
Merge:     112
Counting:  118
Bucket:    109
*/

int main(int argc, char* argv[])
{
    bool outCSV = argc == 2 && std::strcmp(argv[1],"-c") == 0;
    string sortingAlgorithms[] = {"Selection","Heap","Insertion","Bubble","Merge","Counting", "Bucket"};

    if (!outCSV) {
        std::cout<<"Selection Sort, Heap Sort, Insertion Sort, Bubble Sort, Merge Sort, Counting Sort, Bucket Sort\n";
    } else {
        std::cout<<",Time(μs)";
        std::cout<<"Trial,Selection,Heap,Insertion,Bubble,Merge,Counting,Bucket";
    }
    int arr[arraySize];
    for (int i = 1; i <= arraySize; i++) {
        arr[i - 1] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    for (int trial = 0; trial < 10; trial++) {
        if (!outCSV) {
            std::cout << "\n==== TRIAL #" << trial + 1 << " ===="<<std::endl;
        } else {
            std::cout << "\n" << trial + 1;
        }
        for (int i = 0; i < sizeof(sortingAlgorithms)/sizeof(string); i++)
        {
            for (int s = 0; s < 5; s++)
            {
                std::shuffle(std::begin(arr), std::end(arr), g);
            }
            
            runAlgorithm(arr, sortingAlgorithms[i], outCSV);
        }
    }
    std::cout << std::flush;
    return 0;
}

/*
Magic Numbers:
184 - Selection
180 - Quick
184 
174
164
180
177
184
*/