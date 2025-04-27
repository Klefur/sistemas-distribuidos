#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

// Función para generar un número aleatorio en el rango [-1, 1]
double random_double() {
    return (double)rand() / RAND_MAX * 2.0 - 1.0;
}

int main(int argc, char* argv[]) {
    int total_points = 1000000; // Número total de puntos por defecto
    int rank, size;
    int thread_count = 4;

    // Inicialización de MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Verificar argumentos
    if (argc > 1) {
        total_points = atoi(argv[1]);
    }

    if (argc > 2) {
        thread_count = atoi(argv[2]);
    }

    int points_per_process = total_points / size;
    int local_count = 0;

    // Medir el tiempo de ejecución local
    double start_time = MPI_Wtime();

    // Configurar OpenMP
    #pragma omp parallel num_threads(thread_count)
    {
        unsigned int seed = rank * omp_get_thread_num(); // Semilla única por hilo

        #pragma omp for
        for (int i = 0; i < points_per_process; ++i) {
            double x = random_double();
            double y = random_double();
            #pragma omp critical
            {
                if (x * x + y * y <= 1.0) {
                    local_count++;
                }
            }
        }
    }

    double end_time = MPI_Wtime();
    double local_execution_time = end_time - start_time;

    // Reducir resultados usando MPI_Reduce
    int global_count;
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reducir tiempos de inicio y finalización para calcular el tiempo total
    double global_start_time, global_end_time;
    MPI_Reduce(&start_time, &global_start_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&end_time, &global_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // El proceso principal calcula π y muestra el tiempo total
    if (rank == 0) {
        double pi = 4.0 * global_count / total_points;
        double total_execution_time = global_end_time - global_start_time;

        printf("Valor de π: %f\n", pi);
        printf("Tiempo total de ejecución del programa: %f segundos\n", total_execution_time);
    }

    // Finalizar MPI
    MPI_Finalize();

    return 0;
}
