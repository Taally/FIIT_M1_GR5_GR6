#include <iostream>
#include <random>
#include <vector>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <string>

using namespace std;

/*
A: �����-�����������
a00 0   0   0
a10 a11 0   0  
a20 a21 a22 0
a30 a31 a32 a33

B: ������������ (�� ��� �����������������)
b00 b01 b02 b03
0   b11 b12 b13
0   0   b22 b23
0   0   0   b33

�������� �� ��������:
a00 a10 a20 a30 a11 a21 a31 a22 a32 a33
b00 b01 b11 b02 b12 b13 b03 b13 b23 b33
*/

typedef float elem;
//typedef double elem;

const int g_nNumberOfThreads = 8;
const int N = 2880; //2880; // matrix size 
const elem eps = 1e-1;
vector<elem> A; // LowerTriangle matrix A
vector<elem> B; // Simmetric matrix B
vector<elem> realC;

elem random()
{
	elem max = 10;
	return static_cast<elem>(rand() / 1000);
	return static_cast<elem>(rand()) / (static_cast<elem>(RAND_MAX / max));
}

void generateMatrices()
{
	A.resize(N*N);
	B.resize(N*N);

	for (int j = 0; j < N; ++j)
		for (int i = 0; i < N; ++i)
		{
			int index = i + j * N;
			if (j <= i)
			{
				A[index] = random();
				B[index] = random();
			}
			else 
			{
				A[index] = 0;
				B[index] = B[i * N + j];
			}
		}
}

vector<elem> fillTriang()
{
	int nonZeros = (N*N + N) / 2;
	vector<elem> matrix(nonZeros);

	for (int i = 0; i < nonZeros; ++i) {
		
		matrix[i] = random();
	}
	return matrix;
}

bool isZeroLower(int i, int j)
{
	return (i < j);
}

// ������ �������� (i, j) ���������������� ������� � ���������� �������
int getLowerElemNumber(int i, int j, int size)
{
	if (isZeroLower(i, j))
		return -1;
	int k = 0;
	for (int j1 = 0; j1 < j; ++j1)
		k += (size - j1);	
	k += (i - j);
	return k;
}

// ������ �������� (i, j) ������������ ������� � ���������� �������
int getSimElemNumber(int i, int j)
{
	if (i > j)
	{
		int t = i;
		i = j;
		j = t;
	}
	int k = 0;
	for (int j1 = 1; j1 <= j; ++j1)
		k += j1;
	k += i;
	return k;
}

// ������� (i, j) ����� (bI, bJ) � ���������� ������� ���������������� �������
elem elemInLowerBlock(int i, int j, int bI, int bJ, int M)
{
	int iMatr = i + bI * M;
	int jMatr = j + bJ * M;
	if (iMatr < jMatr)
		return 0;
	else
	{
		int idx = getLowerElemNumber(iMatr, jMatr, N);
		return A[idx];
	}
}

// ������� (i, j) ����� (bI, bJ) � ���������� ������� ������������ �������
elem elemInSimBlock(int i, int j, int bI, int bJ, int M)
{
	int iMatr = i + bI * M;
	int jMatr = j + bJ * M;
	int idx;
	idx = getSimElemNumber(iMatr, jMatr);
	return B[idx];
}

vector<elem> SumMatrices(vector<vector<elem>> & matrices)
{
	int length = matrices[0].size();
	vector<elem> res(length);
	for (int i = 0; i < length; ++i)
	{
		elem sum = 0;
		for (int j = 0; j < matrices.size(); ++j)
			sum += matrices[j][i];
		res[i] = sum;
	}
	return res;
}

vector<elem> getBlockForLower(int bI, int bJ, int M)
{
	vector<elem> block(M*M);
	for (int j = 0; j < M; ++j)
		for (int i = 0; i < M; ++i)
		{
			int elI = bI * M + i;
			int elJ = bJ * M + j;
			block[i + j * M] = A[elI + elJ * N];
			//block[i + j * M] = elemInLowerBlock(i, j, bI, bJ, M);
		}
	return block;
}

vector<vector<elem>> CreateABlocks(int M)
{
	int BLOCK_CNT = N / M; // number of blocks on one side
	int size = (BLOCK_CNT * BLOCK_CNT + BLOCK_CNT) / 2;
	
	vector<vector<elem>> aBlocks(size);
	for (int i = 0; i < BLOCK_CNT; ++i)
		for (int j = 0; j <= i; ++j)
		{ 
			int k = getLowerElemNumber(i, j, BLOCK_CNT);
			aBlocks[k] = getBlockForLower(i, j, M);
		}
	return aBlocks;
}

vector<elem> getBlockForSim(int bI, int bJ, int M)
{
	vector<elem> block(M*M);
	for (int j = 0; j < M; ++j)
		for (int i = 0; i < M; ++i)
		{
			int elI = bI * M + i;
			int elJ = bJ * M + j;
			block[i + j * M] = B[elI + elJ * N];
			//block[i + j * M] = elemInSimBlock(i, j, bI, bJ, M);
 		}
	return block;
}

vector<vector<elem>> CreateBBlocks(int M)
{
	int BLOCK_CNT = N / M; // number of blocks on one side
	int size = (BLOCK_CNT * BLOCK_CNT + BLOCK_CNT) / 2;

	vector<vector<elem>> bBlocks(size);
	for (int i = 0; i < BLOCK_CNT; ++i)
		for (int j = i; j < BLOCK_CNT; ++j)
		{
			int k = getSimElemNumber(i, j);
			bBlocks[k] = getBlockForSim(i, j, M);
		}
	return bBlocks;
}

vector<elem> MultBlocks(vector<elem> & blockA, vector<elem> & blockB, int M)
{
	vector<elem> resMatr(M*M);
//#pragma omp parallel for // parallel 2
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			elem res = 0;
			for (int k = 0; k < M; ++k)
			{
				int indA = i + k * M;
				int indB = k + j * M;
				res += blockA[indA] * blockB[indB];
			}
			resMatr[i + j * M] = res;
		}
	}
	return resMatr;
}

vector<elem> MultBlocksTransposedB(vector<elem> & blockA, vector<elem> & blockB, int M)
{
	vector<elem> resMatr(M*M);
//#pragma omp parallel for // parallel 2
	for (int i = 0; i < M; ++i)
	{
		for (int j = 0; j < M; ++j)
		{
			elem res = 0;
			for (int k = 0; k < M; ++k)
			{
				int indA = i + k * M;
				int indB = j + k * M;
				res += blockA[indA] * blockB[indB];
			}
			resMatr[i + j * M] = res;
		}
	}
	return resMatr;
}

// Calculate result block with coords (cI, cJ)
vector<elem> CalcBlockC(int cI, int cJ, 
						vector<vector<elem>> & aBlocks, 
						vector<vector<elem>> & bBlocks, 
						int M)
{
	const int BLOCK_CNT = N / M;
	vector<vector<elem>> blocks(cI + 1);
	
	for (int k = cI; k > cJ; --k)
	{
		int indA = getLowerElemNumber(cI, k, BLOCK_CNT);
		int indB = getSimElemNumber(k, cJ);
		blocks[k] = MultBlocksTransposedB(aBlocks[indA], bBlocks[indB], M);
	}
	for (int k = cJ; k >= 0; --k)
	{
		int indA = getLowerElemNumber(cI, k, BLOCK_CNT);
		int indB = getSimElemNumber(k, cJ);
		blocks[k] = MultBlocks(aBlocks[indA], bBlocks[indB], M);
	}
	return SumMatrices(blocks);
}

// Construct result from blocks
vector<elem> ConstructC(vector<vector<elem>> & blocks, int M)
{
	int BLOCK_CNT = N / M; // number of blocks on one side
	int TOTAL_BLOCK_CNT = BLOCK_CNT * BLOCK_CNT; // total number of blocks
	vector<elem> res(N*N);
	for (int blockNum = 0; blockNum < TOTAL_BLOCK_CNT; ++blockNum)
	{
		for (int i = 0; i < M; ++i)
			for (int j = 0; j < M; ++j)
			{
				int resI = (blockNum % BLOCK_CNT) * M + i;
				int resJ = (blockNum / BLOCK_CNT) * M + j;
				res[resJ * N + resI] = blocks[blockNum][j*M + i];
			}
	}
	return res;
}

// Simple multiplication
vector<elem> multiplyAB() {
	realC.resize(N*N);

//# pragma omp parallel for
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
		{
			elem res = 0;
			for (int k = 0; k < N; ++k) 
				res += A[i + k * N] * B[k + j * N];
			realC[j * N + i] = res;
		}
	return realC;
}

void saveRes(vector<elem> & C)
{
	realC = C;
}

vector<elem> loadRes()
{
	ifstream file;
	file.open("C.txt");
	int a;
	vector<elem> C(N*N);
	int i = 0;
	while (file >> a)
	{
		C[i] = a;
		i++;
	}
	return C;
}

bool checkC(vector<elem> & C)
{
	if (C.size() != realC.size())
		return false;
	for (int i = 0; i < C.size(); ++i)
	{
		if (fabs(C[i] - realC[i]) > eps) {
			cout << C[i] << ", " << realC[i] << endl;
			return false;
		}
	}
	return true;
}

vector<elem> calculateC(int M)
{
	int BLOCK_CNT = N / M; // number of blocks on one side
	int TOTAL_BLOCK_CNT = BLOCK_CNT * BLOCK_CNT; // total number of blocks
	cout << M << endl;

	vector<vector<elem>> aBlocks = CreateABlocks(M);
	vector<vector<elem>> bBlocks = CreateBBlocks(M);

	vector<vector<elem>> blocksC(TOTAL_BLOCK_CNT);

	clock_t begin = clock();
	//#pragma omp parallel for // parallel 1
	for (int blockI = 0; blockI < BLOCK_CNT; ++blockI) {
		for (int blockJ = 0; blockJ < BLOCK_CNT; ++blockJ) {
	
			vector<elem> curBlock(M*M);
			//vector<vector<elem>> blocks;
			//blocks.resize(blockI);

			int blockK;
			for (blockK = blockI; blockK > blockJ; --blockK)
			{
				int indA = getLowerElemNumber(blockI, blockK, BLOCK_CNT);
				int indB = getSimElemNumber(blockK, blockJ);

				//blocks[blockK] = MultBlocksTransposedB(aBlocks[indA], bBlocks[indB], M);

				#pragma omp parallel for // parallel 2
				for (int i = 0; i < M; ++i)
				{
					for (int j = 0; j < M; ++j)
					{
						elem res = 0;
						for (int k = 0; k < M; ++k)
						{
							res += aBlocks[indA][i + k * M] * bBlocks[indB][j + k * M];
						}
						curBlock[i + j * M] += res;
					}
				}
			}
			for (; blockK  >= 0; --blockK)
			{
				int indA = getLowerElemNumber(blockI, blockK, BLOCK_CNT);
				int indB = getSimElemNumber(blockK, blockJ);

				//blocks[blockK] = MultBlocks(aBlocks[indA], bBlocks[indB], M);
				#pragma omp parallel for // parallel 2
				for (int i = 0; i < M; ++i)
				{
					for (int j = 0; j < M; ++j)
					{
						elem res = 0;
						for (int k = 0; k < M; ++k)
						{
							res += aBlocks[indA][i + k * M] * bBlocks[indB][k + j * M];
						}
						curBlock[i + j * M] += res;
					}
				}
			}

			//blocksC[blockJ*BLOCK_CNT + blockI] == SumMatrices(blocks);
			blocksC[blockJ*BLOCK_CNT + blockI] = curBlock;
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / 1000;

	vector<elem> C = ConstructC(blocksC, M);
	bool correct = checkC(C);
	cout << (correct ? "Yes" : "No") << endl;
	cout << elapsed_secs << endl << endl;
	return C;
}

void measureTime()
{
	vector<int> blocksizes = { 1, 6, 10, 15, 20, 24, 30, 36, 40, 60, 72, 80, 96, 120, 144, 160, 180, 240, 360, 480, 720 };
	//vector<int> blocksizes = { 80, 96, 120, 144, 160, 180 };
	generateMatrices();

	clock_t begin = clock();
	multiplyAB();
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / 1000;
	cout << "Simple multiplication took " << elapsed_secs << endl;
	
	cout << "Start calc blocks" << endl;
	int k = blocksizes.size() - 1;
	while (k >= 0) {
		calculateC(blocksizes[k]);
		k = k - 1;
	}
}

int main()
{
	measureTime();
	int wait;
	cin >> wait;
}