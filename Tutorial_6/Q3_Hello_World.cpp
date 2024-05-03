#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>
#include <unistd.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char machine_name[256];
    gethostname(machine_name, 255);

    if (rank == 0) {
        std::cout << "Hello world from master process " << rank << " running on " << machine_name << std::endl;

        std::vector<std::string> messages(size - 1);
        MPI_Status status;
        for (int i = 1; i < size; ++i) {
            MPI_Recv(&messages[i - 1][0], 256, MPI_CHAR, i, 99, MPI_COMM_WORLD, &status);
            std::cout << "Message from process " << i << ": " << messages[i - 1] << std::endl;
        }
    } else {
        char message[256];
        printf(message, "Hello world from process %d running on %s", rank, machine_name);
        MPI_Send(&message[0], 256, MPI_CHAR, 0, 99, MPI_COMM_WORLD); 
    } 

    MPI_Finalize();
    return 0;
}