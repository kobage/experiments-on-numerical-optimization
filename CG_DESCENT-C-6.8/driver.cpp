#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "cg_user.h"
#include "cg_descent.h"
#include "cg_blas.h"
#include <time.h>
#include "testFunctions.h"

void run_cg_descent_Tests(struct testFunction* t, INT number)
{
	makeTestCollection(test);
	clock_t st, fin;
	int repNumber = 10;
	double running_times[10];
	cg_stats Stats;

	puts("The following tests are available. plz enter test number");
	INT i;
	for (i = 0; i < number; i++)
		printf("%d. %s, size of test: %d\n ", i + 1, test[i].name, test[i].size);

	scanf("%d", &i);
	getchar();
	i--;
	if (i < 1 || i > 45) return;

	n = test[i].size;
	printf("Chosen: %s, size of test: %d\n ", test[i].name, test[i].size);

	double *x;
	x = (double *)malloc(n * sizeof(double));
	int j;
	for (j = 0; j < repNumber; j++)
	{
		test[i].StartingGess(x, n);
		cg_parameter Parm;
		cg_default(&Parm);
		//Parm.PrintLevel = 0 ;
		st = clock();

		cg_descent(x, n, &Stats, &Parm, 1.e-6, test[i].value, test[i].grad, test[i].valgrad, NULL);
		//cg_descent (x, n, &Stats, NULL, 1.e-6, test[i].value, test[i].grad, test[i].valgrad, NULL) ;
		fin = clock();

		running_times[j] = ((double)(fin - st)) / CLOCKS_PER_SEC;
	}
	iSort(running_times, 10);

	freopen("results.txt", "w", stdout);
	{
		printf("Median time:   %f\n", running_times[repNumber / 2]);
		printf("Function's value:   %f\n", Stats.f);
		int quarter = n / 4;
		int ii;
		for (ii = 0; ii < quarter; ii++)
		{
			printf("x[%3d] = %-15.11f     x[%3d] = %15.11f     x[%3d] = %15.11f     x[%3d] = %2.11f  \n",
				ii, x[ii], quarter + ii, x[quarter + ii], 2 * quarter + ii, x[2 * quarter + ii], 3 * quarter + ii, x[3 * quarter + ii]);
		}
	}
	free(x);
}

int main(void)
{
	run_cg_descent_Tests(test, TestNumber);
}