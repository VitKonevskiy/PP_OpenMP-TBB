#include <omp.h>
#include <stdio.h>
#include <mpi.h>
#include <Windows.h>
#include <time.h>
#include <mutex>
#define N_THREADS 4
//#define N_PROC 2
#define N 25000000
#define RAND_MAX_DOUBLE 1000

#define PRINT_ARRAY 0
#define PRINT_PARTS_OF_ARRAY_BEFORE_SORT 0
#define PRINT_PARTS_OF_ARRAY_AFTER_SORT 0
#define PRINT_AFTER_COLLECT 0


#define COMPARE_LSD_SORT_WITH_QSORT 1
#define RUN_STEP_BY_STEP 1
#define PRINT_PARALLEL_LSD_SORT_TIME 1
//MPI_Wtime()

// Сортировка подсчетом для типа double по i-му байту
void CountingSort(double* arr_inp, double* arr_out, int size_arr, int byte_num, int start_index, int part_of_thread, int *counter, int offset, int j, int tmp)
{
	// Идея следующая:
	// Каждый байт может иметь 256 различных состояний (диапазон чисел от 0 до 255)
	// Пусть есть байты А и В с номерами состояний s_A и s_B соответственно
	// Байт A больше байта B, если s_A > s_B
	unsigned char* mas = (unsigned char*)arr_inp;

	//int counter[256];// Возможно 256 различных состояний одного байта
	//int offset;
	for (j = 0; j < 256; j++)
	{
		counter[j] = 0;
	}
	//memset(counter, 0, sizeof(int) * 256); // Обнуляем массив
	// mas возвращает символ с кодом от 0 до 255 (иными словами номер состояния байта)
	for (int i = 0; i < part_of_thread; i++)
		counter[mas[8 * i + byte_num]]++; // counter показывает, сколько чисел типа double содержит определенный разряд

										  // Теперь ищем номер состояния байта byte_num, который присутствует в каких-либо элементах double 
										  //int j = 0;
	j = 0;
	for (j; j < 256; j++)
		if (counter[j] != 0)
			break;

	offset = counter[j];// Теперь offset показывает, сколько имеется элементов с определенным байтом (чтобы определить, сколько ячеек массива arr_out уйдет под числа,
						// содержащих байт с номером состояния j)
	counter[j] = 0;// Это характеризует смещение элементов, содержащих байт с номером состояния j. Причем такие элементы будут иметь "наименьший байт" и будут записаны в начале массива arr_out
	j++;

	// Далее считаем смещения и записываем их в counter

	for (j; j < 256; j++)
	{
		tmp = counter[j];
		counter[j] = offset;
		offset += tmp;
	}

	for (int i = 0; i < part_of_thread; i++)
	{
		arr_out[counter[mas[8 * i + byte_num]]] = arr_inp[i];// counter содержит всю необходимую информацию по корректному раскидыванию элементов
		counter[mas[8 * i + byte_num]]++;// Увеличиваем смещение на 1 для элемента (чтобы корректно его записать в ячейку массива arr_out)
	}
}


void LSDSortDouble(double* arr_inp, int start_index, int part_of_thread, int size_arr, double * arr_inp_tmp, double* arr_out_tmp, int *counter, int offset, int j, int tmp)
{
	// Создадим два массива, которые будут содержать по отдельности отрицательные и положительные элементы
	// Для каждого из них по 8 раз запустим соответсвующую сортировку подсчетом
	// После этого сливаем эти массивы и получаем результат


	for (int i = start_index; i < start_index + part_of_thread; i++)
		arr_inp_tmp[i - start_index] = arr_inp[i];

	// Сортируем массив

	CountingSort(arr_inp_tmp, arr_out_tmp, size_arr, 0, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_out_tmp, arr_inp_tmp, size_arr, 1, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_inp_tmp, arr_out_tmp, size_arr, 2, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_out_tmp, arr_inp_tmp, size_arr, 3, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_inp_tmp, arr_out_tmp, size_arr, 4, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_out_tmp, arr_inp_tmp, size_arr, 5, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_inp_tmp, arr_out_tmp, size_arr, 6, start_index, part_of_thread, counter, offset, j, tmp);
	CountingSort(arr_out_tmp, arr_inp_tmp, size_arr, 7, start_index, part_of_thread, counter, offset, j, tmp);
	//Записываем результат
	for (int i = start_index; i < start_index + part_of_thread; i++)
		arr_inp[i] = arr_inp_tmp[i - start_index];

}


// Сортировка подсчетом для типа double по i-му байту
void CountingSortStepByStep(double* arr_inp, double* arr_out, int size_arr, int byte_num)
{
	// Идея следующая:
	// Каждый байт может иметь 256 различных состояний (диапазон чисел от 0 до 255)
	// Пусть есть байты А и В с номерами состояний s_A и s_B соответственно
	// Байт A больше байта B, если s_A > s_B
	unsigned char* mas = (unsigned char*)arr_inp;

	int counter[256];// Возможно 256 различных состояний одного байта
	int offset;

	memset(counter, 0, sizeof(int) * 256); // Обнуляем массив

										   // mas возвращает символ с кодом от 0 до 255 (иными словами номер состояния байта)
	for (int i = 0; i < size_arr; i++)
		counter[mas[8 * i + byte_num]]++; // counter показывает, сколько чисел типа double содержит определенный разряд

										  // Теперь ищем номер состояния байта byte_num, который присутствует в каких-либо элементах double 
	int j = 0;
	for (j; j < 256; j++)
		if (counter[j] != 0)
			break;

	offset = counter[j];// Теперь offset показывает, сколько имеется элементов с определенным байтом (чтобы определить, сколько ячеек массива arr_out уйдет под числа,
						// содержащих байт с номером состояния j)
	counter[j] = 0;// Это характеризует смещение элементов, содержащих байт с номером состояния j. Причем такие элементы будут иметь "наименьший байт" и будут записаны в начале массива arr_out
	j++;

	// Далее считаем смещения и записываем их в counter

	for (j; j < 256; j++)
	{
		int tmp = counter[j];
		counter[j] = offset;
		offset += tmp;
	}

	for (int i = 0; i < size_arr; i++)
	{
		arr_out[counter[mas[8 * i + byte_num]]] = arr_inp[i];// counter содержит всю необходимую информацию по корректному раскидыванию элементов
		counter[mas[8 * i + byte_num]]++;// Увеличиваем смещение на 1 для элемента (чтобы корректно его записать в ячейку массива arr_out)
	}
}


void LSDSortDoubleStepByStep(double* arr_inp, int size_arr)
{
	// Создадим два массива, которые будут содержать по отдельности отрицательные и положительные элементы
	// Для каждого из них по 8 раз запустим соответсвующую сортировку подсчетом
	// После этого сливаем эти массивы и получаем результат

	double* arr_inp_tmp;
	double* arr_out_tmp;
	int size_arr_plus = 0,
		size_arr_minus = 0;

	arr_inp_tmp = new double[size_arr];
	arr_out_tmp = new double[size_arr];

	for (int i = 0; i < size_arr; i++)
		arr_inp_tmp[i] = arr_inp[i];

	// Сортируем массив
	CountingSortStepByStep(arr_inp_tmp, arr_out_tmp, size_arr, 0);
	CountingSortStepByStep(arr_out_tmp, arr_inp_tmp, size_arr, 1);
	CountingSortStepByStep(arr_inp_tmp, arr_out_tmp, size_arr, 2);
	CountingSortStepByStep(arr_out_tmp, arr_inp_tmp, size_arr, 3);
	CountingSortStepByStep(arr_inp_tmp, arr_out_tmp, size_arr, 4);
	CountingSortStepByStep(arr_out_tmp, arr_inp_tmp, size_arr, 5);
	CountingSortStepByStep(arr_inp_tmp, arr_out_tmp, size_arr, 6);
	CountingSortStepByStep(arr_out_tmp, arr_inp_tmp, size_arr, 7);

	//Записываем результат
	for (int i = 0; i < size_arr; i++)
		arr_inp[i] = arr_inp_tmp[i];
	
	delete[]arr_inp_tmp;
	delete[]arr_out_tmp;
}

void CollectThreads(double *array, double *tmpArray, int start_index, int part_of_thread, int curr_step, int a, int b, int Size1, int Size2, int j)
{
	a = start_index;
	b = start_index + part_of_thread;
	Size1 = start_index + part_of_thread;
	Size2 = b + part_of_thread;
	j = start_index;
	if ((omp_get_thread_num() % curr_step) == 0)				//При выборе curr step =2 :  0  2  4  6  8 
	{
		while ((a != Size1) && (b != Size2))
		{
			if (array[a] <= array[b])
			{
				tmpArray[j] = array[a];
				a++;
			}
			else
			{
				tmpArray[j] = array[b];
				b++;
			}
			j++;
		}
		//Дописываем остаток
		if (a == Size1)
			for (j = b; j < Size2; j++)
				tmpArray[j] = array[j];
		else
			for (j = a; j < Size1; j++)
				tmpArray[j + part_of_thread] = array[j];
	}
}

/*double *CollectProc(double *rbuf,int partOfProc,int curr_step,int rank)
{
	double *array ;
	double *tmpArray;
	if (rank %curr_step == 0)
	{
		tmpArray = new double[partOfProc];
		MPI_Recv(tmpArray, partOfProc, MPI_DOUBLE, rank + (curr_step / 2), 0, MPI_COMM_WORLD, 0);
		array = new double[partOfProc* curr_step];
		

		int a = 0;
		int b = 0;
		int j = 0;

		while ((a != partOfProc) && (b != partOfProc))
		{
			if (rbuf[a] <= tmpArray[b])
			{
				array[j] = rbuf[a];
				a++;
			}
			else
			{
				array[j] = tmpArray[b];
				b++;
			}
			j++;
		}
		//Дописываем остаток
		if (a == partOfProc)
			for (j = b; j < partOfProc; j++)
				array[a + j] = tmpArray[j];
		else
			for (j = a; j < partOfProc; j++)
				array[j + b] = rbuf[j];

		if (rank %curr_step == 0)
		{
			delete[]tmpArray;
		}

		return array;
	}
	else
	{
		MPI_Send(rbuf,partOfProc,MPI_DOUBLE,rank - (curr_step/2),0,MPI_COMM_WORLD);
		return NULL;
	}

	
}
*/

double* CollectProc(double *array0, int size0, double *array1, int size1)
{

	double *result = new double[size0 + size1];
	int i = 0, j = 0;
	int index = 0;

	while ((i<size0) && (j<size1))
	{
		if (array0[i] < array1[j])
		{
			result[index] = array0[i];
			i++;
		}
		else
		{
			result[index] = array1[j];
			j++;
		}

		index++;
	}

	while (i < size0)
	{
		result[index] = array0[i];
		index++;
		i++;
	}

	while (j < size1)
	{
		result[index] = array1[j];
		index++;
		j++;
	}
	delete[]array0;
	delete[]array1;
	return result;
}

static int compare(const void * a, const void * b)
{
	if (*(double*)a > *(double*)b) return 1;
	else if (*(double*)a < *(double*)b) return -1;
	else return 0;
}


int main(int argc,char *argv[])
{
	srand(time(0));

	int rank=0;
	int numtasks = 0;
	int root = 0;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	int N_PROC = numtasks;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int partOfProc = N / numtasks;
	/*if (rank == 0)
	{
		printf("rank = %d", rank);
		fflush(stdout);
	}*/
	int count=0;
/*#pragma omp parallel num_threads(N_THREADS)
	{
		count++;
		/*if ((rank == 0) && (omp_get_thread_num()==0))
		{
			printf("r %d, t %d\n", rank, omp_get_thread_num());
			fflush(stdout);
		}
	}
	if (rank == 0)
	{
		count--;
		printf("count = %d", count);
	}*/
	
	double *array =new double[N];
	double *arrayForStepByStep = new double[N];
	double *rbuf = new double[partOfProc];
	double *tmpArray = new double[partOfProc];
	double startParallelTime=0, endParallelTime = 0;

	if (rank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			array[i] = (double)(rand()) / RAND_MAX * RAND_MAX_DOUBLE;
			arrayForStepByStep[i] = array[i];
			if (PRINT_ARRAY)
			{
				printf("%lf \n", array[i]);
				fflush(stdout);
			}
		}
	}


	//--------------------------------------------- START OF PARALLEL PART OF PROGRAMM------------------------------------------------//
	if (rank == 0)
	{
		startParallelTime = MPI_Wtime();
	}
	//Раздаем всем по кусочку
	MPI_Scatter(array, partOfProc, MPI_DOUBLE, rbuf, partOfProc, MPI_DOUBLE, root, MPI_COMM_WORLD);
	
	if (PRINT_PARTS_OF_ARRAY_BEFORE_SORT)
	{
		for (int i = 0; i < partOfProc; i++)
		{
			printf("BeforeSort! rank = %d, rbuf[%d] = %lf\n", rank,i, rbuf[i]);
			fflush(stdout);
		}
	}

	//------------------------------------------- START OPEN MP-------------------------------------------//
#pragma omp parallel num_threads(N_THREADS)
	{
		int offset = 0;
		int *counter = new int[256];
		int partOfThread = partOfProc / omp_get_num_threads();
		int start_index = partOfThread * omp_get_thread_num();
		int numthread = omp_get_thread_num();
		int j = 0;
		int tmp = 0;
		double *arr_out_tmp = new double[partOfThread];
		double *arr_inp_tmp = new double[partOfThread];

		//переменные необходимые для слияния
		int a = start_index;
		int b = 0;
		int Size1 = a + partOfThread;
		int Size2 = b + partOfThread;


		LSDSortDouble(rbuf, start_index, partOfThread, partOfProc, arr_inp_tmp, arr_out_tmp, counter, offset, j, tmp);

#pragma omp barrier
		CollectThreads(rbuf, tmpArray, start_index, partOfThread, 2, a, b, Size1, Size2, j);			//После выполнения функции, ведущие потоки :   0	2	4	6

#pragma omp barrier
		CollectThreads(tmpArray, rbuf, start_index, 2 * partOfThread, 4, a, b, Size1, Size2, j);		//После выполнения функции, ведущие потоки :   0	4	
//#pragma omp barrier
		
		delete[]arr_out_tmp;
		delete[]arr_inp_tmp;
		delete[]counter;
	}
	//--------------------------------------------- END OPEN MP----------------------------------------------//

	if (PRINT_PARTS_OF_ARRAY_AFTER_SORT)
	{
		for (int i = 0; i < partOfProc; i++)
		{
			printf("AfterSort! rank = %d, rbuf[%d] = %lf\n", rank, i, rbuf[i]);
			fflush(stdout);
		}
	}

	int step = numtasks;
	
	MPI_Status status;
	int m = 1;
	int arr0_size = partOfProc;
	//-----------------------------START COLLECT IN PROC------------------------------//
	while (step > 1)
	{
		step = step / 2;
		if ((rank  - m) % (2 * m) == 0)
		{

			MPI_Send(&arr0_size, 1, MPI_INT, rank - m, 0, MPI_COMM_WORLD);
			MPI_Send(rbuf, arr0_size, MPI_DOUBLE, rank - m, 0, MPI_COMM_WORLD);
			//printf("SEND from %d, TO %d", rank, rank - m);
			fflush(stdout);
		}
		if ((rank % (2 * m) == 0) && (numtasks - rank > m))
		{
			int arr1_size;
			double *tmpArray;
			//printf("RECV from %d, IN ME %d", rank + m, rank);
			fflush(stdout);
			MPI_Recv(&arr1_size, 1, MPI_INT, rank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tmpArray = new double[arr1_size];
			MPI_Recv(tmpArray, arr1_size, MPI_DOUBLE, rank + m, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			rbuf = CollectProc(rbuf, arr0_size, tmpArray, arr1_size);			//Память от передаваемого в CollectProc rbuf и tmpArray удаляется внутри функции.

			arr0_size = arr0_size + arr1_size;
		}
		m = 2 * m;
	}
	//-----------------------------END COLLECT IN PROC--------------------------------//
	

	//----------------------------------------------------------- END OF PARALLEL PART OF PROGRAMM-------------------------------------------------------------//
	if (rank == 0)
	{
		endParallelTime = MPI_Wtime();
	}
	if (rank == 0)
	{
		if (RUN_STEP_BY_STEP)
		{
			double startTime = MPI_Wtime();
			LSDSortDoubleStepByStep(arrayForStepByStep, N);
			double endTime = MPI_Wtime();
			printf("time of StepByStepLSD_Sort_Double = %lf\n", endTime - startTime);
		}
		if (PRINT_AFTER_COLLECT)
		{
			for (int i = 0; i < N; i++)
			{
				printf("result[%d] = %lf\n",i,rbuf[i]);
				fflush(stdout);
			}
		}
		if (COMPARE_LSD_SORT_WITH_QSORT)
		{
			bool flag1 = true;
			bool flag2 = true;
			double startQsort = MPI_Wtime();
			qsort(array, N, sizeof(double), compare);
			double endQsort = MPI_Wtime();
			for (int i = 0; i < N; i++)
			{
				if (rbuf[i] != array[i])
					flag1 = false;
				if (RUN_STEP_BY_STEP)
				{
					if (arrayForStepByStep[i] != array[i])
						flag2 = false;
				}
			}
			printf("Proc + thread LSD Sort ? std::qsort  %d\n", flag1);
			
			if (RUN_STEP_BY_STEP)
			{
				printf("   StepByStep LSD Sort ? std::qsort  %d\n", flag2);
			}

			printf("time of std::qsort =%lf \n", endQsort - startQsort);

		}
		if (PRINT_PARALLEL_LSD_SORT_TIME)
		{
			printf("time of ParallelLSD_Sort_Double = %lf\n", endParallelTime - startParallelTime);
		}

	}


	
	MPI_Barrier(MPI_COMM_WORLD);
	delete[]tmpArray;
	delete[]array;
	delete[]rbuf;
	delete[]arrayForStepByStep;
	
	MPI_Finalize();
	return 0;
}