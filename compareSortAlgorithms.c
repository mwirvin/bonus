#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

int extraMemoryAllocated;

void *Alloc(size_t sz)
{
	extraMemoryAllocated += sz;
	size_t* ret = malloc(sizeof(size_t) + sz);
	*ret = sz;
	printf("Extra memory allocated, size: %ld\n", sz);
	return &ret[1];
}

void DeAlloc(void* ptr)
{
	size_t* pSz = (size_t*)ptr - 1;
	extraMemoryAllocated -= *pSz;
	printf("Extra memory deallocated, size: %ld\n", *pSz);
	free((size_t*)ptr - 1);
}

size_t Size(void* ptr)
{
	return ((size_t*)ptr)[-1];
}
void merge(int pData[], int l, int m, int r){
    int i, j, k;
    int n1 = m - l + 1; // num of elements in first subarry
    int n2 = r - m;     // second array

    // dma for arrays
    int* L = (int*)Alloc(n1 * sizeof(int));
    int* R = (int*)Alloc(n2 * sizeof(int));
    // copy data
    for (i = 0; i < n1; i++)
        L[i] = pData[l + i];
    for (j = 0; j < n2; j++)
        R[j] = pData[m + 1+ j];


	// needed inital value
    i = 0; 
    j = 0; 
    k = l;  // merged array
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            pData[k] = L[i];
            i++;
        } else {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    // copy elements of L array
    while (i < n1) {
        pData[k] = L[i];
        i++;
        k++;
    }

    // copy elements of R array
    while (j < n2) {
        pData[k] = R[j];
        j++;
        k++;
    }

    DeAlloc(L);
    DeAlloc(R);
}
void heapify(int arr[], int n, int i) {
    int largest = i; // Initialize largest as root
    int left = 2*i + 1; // left = 2*i + 1
    int right = 2*i + 2; // right = 2*i + 2

    // If left child is larger than root
    if (left < n && arr[left] > arr[largest])
        largest = left;

    // If right child is larger than largest so far
    if (right < n && arr[right] > arr[largest])
        largest = right;

    // If largest is not root
    if (largest != i) {
        int swap = arr[i];
        arr[i] = arr[largest];
        arr[largest] = swap;

        // Recursively heapify the affected sub-tree
        heapify(arr, n, largest);
    }
}

// implements heap sort
// extraMemoryAllocated counts bytes of memory allocated
void heapSort(int arr[], int n) {
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(arr, n, i);

    // One by one extract an element from heap
    for (int i=n-1; i>0; i--) {
        // Move current root to end
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        // call max heapify on the reduced heap
        heapify(arr, i, 0);
    }
}
// implement merge sort
// extraMemoryAllocated counts bytes of extra memory allocated

void merge(int pData[], int l, int m, int r){
    int i, j, k;
    int n1 = m - l + 1; // num of elements in first subarry
    int n2 = r - m;     // second array

    // dma for arrays
    int* L = (int*)Alloc(n1 * sizeof(int));
    int* R = (int*)Alloc(n2 * sizeof(int));
    // copy data
    for (i = 0; i < n1; i++)
        L[i] = pData[l + i];
    for (j = 0; j < n2; j++)
        R[j] = pData[m + 1+ j];


	// needed inital value
    i = 0; 
    j = 0; 
    k = l;  // merged array
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            pData[k] = L[i];
            i++;
        } else {
            pData[k] = R[j];
            j++;
        }
        k++;
    }

    // copy elements of L array
    while (i < n1) {
        pData[k] = L[i];
        i++;
        k++;
    }

    // copy elements of R array
    while (j < n2) {
        pData[k] = R[j];
        j++;
        k++;
    }

    DeAlloc(L);
    DeAlloc(R);
}

void mergeSort(int pData[], int l, int r)
{
    if (l < r) {
        int m = l+(r-l)/2; // new mid point

        
        mergeSort(pData, l, m);
        mergeSort(pData, m+1, r);

		// merge sorted halves
        merge(pData, l, m, r);
    }
}
// implement insertion sort
// extraMemoryAllocated counts bytes of memory allocated
void insertionSort(int* pData, int n) {
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = pData[i]; // The element to be positioned
        j = i - 1;

        // Move elements of pData[0..i-1], that are greater than key,
        // to one position ahead of their current position
        while (j >= 0 && pData[j] > key) {
            pData[j + 1] = pData[j];
            j = j - 1;
        }
        pData[j + 1] = key;
    }
}
void bubbleSort(int* pData, int n)
{
    int i, j, temp;
    for (i = 0; i < n-1; i++)     
    {
        // Last i elements are already in place  
        for (j = 0; j < n-i-1; j++) 
        {
            if (pData[j] > pData[j+1])
            {
                temp = pData[j];
                pData[j] = pData[j+1];
                pData[j+1] = temp;
            }
        }
    }
}
void selectionSort(int* pData, int n)
{
    int i, j, minIdx, temp;

    for (i = 0; i < n-1; i++)
    {
        minIdx = i;
        for (j = i+1; j < n; j++)
        {
            if (pData[j] < pData[minIdx])
            {
                minIdx = j;
            }
        }
        
        // Swap the found minimum element with the first element
        temp = pData[minIdx];
        pData[minIdx] = pData[i];
        pData[i] = temp;
    }
}


// parses input file to an integer array
int parseData(char *inputFileName, int **ppData)
{
	FILE* inFile = fopen(inputFileName,"r");
	int dataSz = 0;
	int i, n, *data;
	*ppData = NULL;
	
	if (inFile)
	{
		fscanf(inFile,"%d\n",&dataSz);
		*ppData = (int *)malloc(sizeof(int) * dataSz);
		// Implement parse data block
		if (*ppData == NULL)
		{
			printf("Cannot allocate memory\n");
			exit(-1);
		}
		for (i=0;i<dataSz;++i)
		{
			fscanf(inFile, "%d ",&n);
			data = *ppData + i;
			*data = n;
		}

		fclose(inFile);
	}
	
	return dataSz;
}
// prints first and last 100 items in the data array
void printArray(int pData[], int dataSz)
{
	int i, sz = dataSz - 100;
	printf("\tData:\n\t");
	for (i=0;i<100;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\t");
	
	for (i=sz;i<dataSz;++i)
	{
		printf("%d ",pData[i]);
	}
	printf("\n\n");
}

int main(void)
{
	clock_t start, end;
	int i;
    double cpu_time_used;
	// char* fileNames[] = {"input1.txt", "input2.txt", "input3.txt"};
	char* fileNames[] = {"input1.txt"};
	
	
	for (i=0;i<3;++i)
	{
		int *pDataSrc, *pDataCopy;
		int dataSz = parseData(fileNames[i], &pDataSrc);
		
		if (dataSz <= 0)
			continue;
		
		pDataCopy = (int *)Alloc(sizeof(int)*dataSz);
	
		printf("---------------------------\n");
		printf("Dataset Size : %d\n",dataSz);
		printf("---------------------------\n");
		
		printf("Selection Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		selectionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Insertion Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		insertionSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

		printf("Bubble Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		bubbleSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		printf("Merge Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		mergeSort(pDataCopy, 0, dataSz - 1);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);

                printf("Heap Sort:\n");
		memcpy(pDataCopy, pDataSrc, dataSz*sizeof(int));
		extraMemoryAllocated = 0;
		start = clock();
		heapSort(pDataCopy, dataSz);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("\truntime\t\t\t: %.1lf\n",cpu_time_used);
		printf("\textra memory allocated\t: %d\n",extraMemoryAllocated);
		printArray(pDataCopy, dataSz);
		
		DeAlloc(pDataCopy);
		DeAlloc(pDataSrc);
	}
	
}
