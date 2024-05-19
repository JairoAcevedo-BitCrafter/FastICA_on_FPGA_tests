#include <array_utilities.h>
#include <matrix.h>
#include <matrix_multiplier_arch.h>
#include <time.h>

#define MAX_16BIT 65535
#define MAX_MACS 128
#define MAX_MAC_SIZE 1024
#define MIN_RANGE -18
#define MAX_RANGE 40
#define MATRIX_TEST 10
#define ERROR_FILE "error.csv"
#define LATENCY_FILE "latency.csv"
#define RELOJ CLOCK_REALTIME

void random_Matrix_generate(float_matrix* matrix,uint32_t rows, uint32_t cols, double max_range){
    for (uint32_t i=0;i<rows;i++){
        for(uint32_t j=0;j<cols;j++){
            float random_number = ((float)random_number_generate( (MAX_16BIT) )) * (float)pow(2.0,max_range - 16.0);
            *(matrix->address+i*matrix->cols+j)=random_number;
        }
    }
}

void random_matrix_size_generate(uint32_t max_macs, uint32_t max_mac_size, int32_t* A_rows, int32_t* comm_dim, int32_t* B_cols){
    *A_rows = random_number_generate(MAX_MACS);
    *A_rows = *A_rows <= 0 ? ((*A_rows * -1) < 2 ? 2 : *A_rows * -1): *A_rows;
    uint32_t max_B_cols = (uint32_t) ((float)MAX_MACS / (float)*A_rows);
    *B_cols = random_number_generate(max_B_cols);
    *B_cols = *B_cols <= 0 ? ((*B_cols * -1) < 1 ? 1 : *B_cols * -1): *B_cols;
    *comm_dim = (int32_t) random_number_generate(MAX_MAC_SIZE);
    *comm_dim = *comm_dim <= 0 ? ((*comm_dim * -1) < 2 ? 2 : *comm_dim * -1): *comm_dim;
}
int32_t random_range(uint32_t limit, uint32_t offset){
    int32_t max_range=abs(random_number_generate(limit))+offset;
    return max_range;
}
void error_calculate_matrix(float_matrix* software_matrix, float_matrix* hardware_matrix, float_matrix* error_matrix){
    //error_matrix[3] = max_error, min_error, prom_error
    //posicion de la matriz resultado es y asi poder luego saber de que producto se origina
    float error= 100 * fabs((*(software_matrix->address) - *(hardware_matrix->address))/ (*(software_matrix->address) == 0.0 ? 1.0 : *(software_matrix->address)));
    *error_matrix->address=error;//max_error
    *(error_matrix->address +1)=error;//min_error
    *(error_matrix->address +2)=0.0;//prom_error
    for (int i=0;i<software_matrix->rows;i++){
        for (int j=0;j<software_matrix->cols;j++){
            error = 100 * fabs(( *(software_matrix->address + i*software_matrix->cols + j) - *(hardware_matrix->address + i*software_matrix->cols + j) ) / (*(software_matrix->address + i*software_matrix->cols + j)==0.0 ? 1.0 : *(software_matrix->address + i*software_matrix->cols + j)) );
            if(error > *(error_matrix->address)){
                *error_matrix->address=error;
            }
            else if(error < *(error_matrix->address +1)){
                *(error_matrix->address +1)=error;
            }
            *(error_matrix->address + 2) += (error/(software_matrix->rows*software_matrix->cols));
        }
    }
}

void write_to_csv_error(char* string, float* error_array, uint32_t A_rows, uint32_t comm_dimm, uint32_t B_cols) {
    // Open the CSV file in append mode
    FILE *file = fopen(string, "a");
    if (file == NULL) {
        // File doesn't exist, create it
        file = fopen(string, "w");
        if (file == NULL) {
            printf("Error creating error file.\n");
            return;
        }
        // Write header line if file is newly created
        fprintf(file, "max_error,min_error,prom_error,rangeA,rangeB\n");
    }
    
    // Write the values to the file
    fprintf(file, "%.6f,%.6f,%.6f,%d,%d,%d\n", error_array[0],error_array[1],error_array[2], A_rows,comm_dimm,B_cols);
    
    // Close the file
    fclose(file);
}
void write_to_csv_latency(char* string, double* latency_array, uint32_t A_rows, uint32_t comm_dimm, uint32_t B_cols) {
    // Open the CSV file in append mode
    FILE *file = fopen(string, "a");
    if (file == NULL) {
        // File doesn't exist, create it
        file = fopen(string, "w");
        if (file == NULL) {
            printf("Error creating latency file.\n");
            return;
        }
        // Write header line if file is newly created
        fprintf(file, "max_latency_SW, min_latency_SW, prom_latency_SW, max_latency_HW, min_latency_HW, prom_latency_HW\n");
    }
    
    // Write the values to the file
    fprintf(file, "%.f,%.f,%.f,%.f,%.f,%.f,%.f,%.f,%.f,%i,%i,%i\n", latency_array[0],latency_array[1],latency_array[2], latency_array[3],latency_array[4],latency_array[5],latency_array[6],latency_array[7],latency_array[8],A_rows,comm_dimm,B_cols);
    
    // Close the file
    fclose(file);
}
float signnum(float value){
    return (value<0? -1 : 1);
}
int vector_abs_comparison(float_matrix *v1, float_matrix *v2, float tol){
    float error;
    if(v1->rows != v2->rows || v1->cols!= v2->cols){
        return 0;
    }
    for(uint16_t i=0; i<v1->cols; i++){
        error=fabs(fabs(*(v1->address+i))-fabs(*(v2->address+i)));
        if(error>tol){
            return 0;
        }
    }
    return 1;
}
int mean_matrix_comparison(float_matrix *mean_matrix, float_matrix *matrix1, float_matrix *matrix2, float tol, float number_unmixing_matrices){
    float_matrix vector1, vector2;
    float_matrix_create(&vector1,1,matrix1->cols);
    float_matrix_create(&vector2,1,matrix2->cols);
    int flag;
    if(matrix1->cols != matrix2->cols || matrix1->rows!= matrix2->rows){
        printf("Size of Matrices arent equal\n");
        return 0;
    }
    for (uint16_t i=0;i<matrix1->rows;i++){
        for (uint16_t j=0;j<matrix2->rows;j++){
            vector1.address=matrix1->address+i*matrix1->cols;
            vector2.address=matrix2->address+j*matrix2->cols;
            flag=vector_abs_comparison(&vector1,&vector2,tol);
            if(flag){
                //mean of vectors
                for(uint16_t k=0; k<matrix1->cols; k++){
                    *(mean_matrix+mean_matrix->rows*i+k)=(signnum(*(vector1.address+k))*fabs(*(vector1.address+k))+fabs(*(vector2.address+k)))/number_unmixing_matrices;
                }
                break;
            }
        }
        if(!flag){
            printf("Matrices are not similar\n");
            return 0;
        }
    }
    return 1;
    float_matrix_destroy(&vector1);
    float_matrix_destroy(&vector2);
}
void mean_unmixing_calculate(float tol, float_matrix *mean_unmixing_matrix, float_matrix *array_unmixing_matrix, uint16_t number_unmixing_matrices){
    for(uint16_t i=0; i<(number_unmixing_matrices-1); i++){
      mean_matrix_comparison(mean_unmixing_matrix,array_unmixing_matrix+i,array_unmixing_matrix+i+1,tol,number_unmixing_matrices); 
    }
}

int main() {
    srand(time(NULL));

    int memory_file;
    //CREATE MULTIPLIER
    matrix_multiplier multiplier1;
    matrix_multiplier_create(&multiplier1, 0, 7, 10, memory_file);

    struct timespec start_SW, end_SW, start_HW_CARGAR, end_HW_CARGAR,start_HW_WAIT, end_HW_WAIT, start_total, end_total;
    struct timespec start_HW_LEER, end_HW_LEER;
    uint32_t iter=0;
    for (int j=1;j<=128;j+=128){
        int p=5;
        while(p<=128/j){
            int i=1;
            while(i<=MAX_MAC_SIZE){
                float error_array[3]; //max_error, min_error, prom_error
                error_array[2]=0.0;//prom_error
                double latency[8]; //max_latency_SW, min_latency_SW, prom_latency_SW, max_latency_HW, min_latency_HW, prom_latency_HW parte 1
                //prom_latency parte 2, prom_latency parte 3
                latency[2]=0.0;
                latency[5]=0.0;
                latency[6]=0.0;
                latency[7]=0.0;
                latency[8]=0.0;
                uint32_t A_rows;
                uint32_t B_cols;
                uint32_t comm_dim;
                A_rows=p;
                B_cols=j;
                comm_dim=i;
                printf("A_rows, comm_dim, B_cols, total es %d %d %d %d\n",A_rows,comm_dim, B_cols, A_rows*B_cols);
                printf("iter es %d\n",iter);
                iter++;
                for(int k=0;k<MATRIX_TEST;k++){
                    int32_t rangeA;
                    int32_t rangeB;
                    rangeA = random_range(MAX_RANGE-MIN_RANGE,MIN_RANGE); 
                    rangeB = random_range(MAX_RANGE-MIN_RANGE,MIN_RANGE); 

                    float_matrix matrixA;
                    float_matrix_create(&matrixA,A_rows,comm_dim);
                    matrixA.max_exponent = rangeA;
                    random_Matrix_generate(&matrixA, A_rows,comm_dim, rangeA);
                    float_matrix matrixB;
                    float_matrix_create(&matrixB,comm_dim, B_cols);
                    matrixB.max_exponent = rangeB;
                    random_Matrix_generate(&matrixB, comm_dim, B_cols, rangeB);
                    
                    double latency_HW_CARGAR;
                    double latency_HW_WAIT;
                    double latency_HW_LEER;
                    //Producto de matrices SW
                    float_matrix software_result;
                    float_matrix_create(&software_result,A_rows,B_cols);
                    clock_gettime(RELOJ, &start_SW);
                    float_matrix_multiply(&matrixA,&matrixB,&software_result);
                    clock_gettime(RELOJ, &end_SW);
                    // Calculate elapsed time in clock ticks
                    double latency_SW;
                    latency_SW = (double) ((double)(end_SW.tv_sec - start_SW.tv_sec)*1000000000.0 + (double)(end_SW.tv_nsec - start_SW.tv_nsec));
                    //Producto de matrices HW
                    double latency_HW_total;
                    float_matrix hardware_result;
                    float_matrix_create(&hardware_result, A_rows, B_cols);
                    clock_gettime(RELOJ, &start_HW_CARGAR);
                    clock_gettime(RELOJ, &start_total);
                    matrix_multiplier_multiply(&matrixA, &matrixB, &multiplier1);
                    clock_gettime(RELOJ, &end_HW_CARGAR);
                    clock_gettime(RELOJ, &start_HW_WAIT);
                    while(!matrix_multiplier_multiply_done(&multiplier1));
                    clock_gettime(RELOJ, &end_HW_WAIT);
                    clock_gettime(RELOJ, &start_HW_LEER);
                    matrix_multiplier_get_result(&hardware_result, &multiplier1);
                    clock_gettime(RELOJ, &end_total);
                    clock_gettime(RELOJ, &end_HW_LEER);
                    latency_HW_total = (double) ((end_total.tv_sec - start_total.tv_sec)*1000000000.0 + (end_total.tv_nsec - start_total.tv_nsec));
                    latency_HW_CARGAR = (double) ((end_HW_CARGAR.tv_sec - start_HW_CARGAR.tv_sec)*1000000000.0 + (end_HW_CARGAR.tv_nsec - start_HW_CARGAR.tv_nsec));
                    latency_HW_WAIT = (double) ((end_HW_WAIT.tv_sec - start_HW_WAIT.tv_sec)*1000000000.0 + (end_HW_WAIT.tv_nsec - start_HW_WAIT.tv_nsec));
                    latency_HW_LEER = (double) ((end_HW_LEER.tv_sec - start_HW_LEER.tv_sec)*1000000000.0 + (end_HW_LEER.tv_nsec - start_HW_LEER.tv_nsec));
                    if(k==0){
                        latency[5]=latency_HW_CARGAR/MATRIX_TEST;
                        latency[6]=latency_HW_WAIT/MATRIX_TEST;
                        latency[7]=latency_HW_LEER/MATRIX_TEST;
                        latency[8]=latency_HW_total/MATRIX_TEST;
                    }
                    else{
                        latency[5] += latency_HW_CARGAR/MATRIX_TEST;
                        latency[6] += latency_HW_WAIT/MATRIX_TEST; 
                        latency[7] += latency_HW_LEER/MATRIX_TEST;
                        latency[8] += latency_HW_total/MATRIX_TEST; 
                    }
                    //printf("Elapsed es %ld y %ld\n", end_HW, start_HW);

                    //COMPROBACIÓN latencias
                    if(k==0){
                        for(int m=0;m<2;m++){
                            latency[m]=latency_SW;
                        }//max_latency_SW, min_latency_SW, prom_latency_SW, max_latency_HW, min_latency_HW, prom_latency_HW
                        for(int n=3;n<=4;n++){
                            latency[n]=latency_HW_total;
                        }
                    }
                    else{
                        if(latency[0] < latency_SW)
                            latency[0]=latency_SW;
                        else if(latency[1] > latency_SW)
                            latency[1]=latency_SW;        
                        if(latency[3]<latency_HW_total)
                            latency[3]=latency_HW_total;
                        else if(latency[4]>latency_HW_total)
                            latency[4]=latency_HW_total;
                    }
                    latency[2] += (latency_SW/MATRIX_TEST);
                    float_matrix_destroy(&matrixA);
                    float_matrix_destroy(&matrixB);
                    float_matrix error_matrix;
                    float_matrix_create(&error_matrix,1,3);
                    error_calculate_matrix(&software_result,&hardware_result,&error_matrix);
                    //float_matrix_print(&error_matrix,' ');
                    if(k==0){
                        error_array[0]=*error_matrix.address;
                        error_array[1]=*(error_matrix.address+1);   
                    }
                    else{
                        if(error_array[0]<*error_matrix.address){
                            error_array[0]=*error_matrix.address;
                        }
                        else if(error_array[1] > *(error_matrix.address+1)){
                            error_array[1]=*(error_matrix.address+1);
                        }
                    }
                    error_array[2] += (*(error_matrix.address+2)/MATRIX_TEST);
                    //float_matrix_print(&software_result, ' ');
                    //float_matrix_print(&hardware_result, ' ');
                    //printf("k es %d\n",k);
                    float_matrix_destroy(&software_result);
                    float_matrix_destroy(&hardware_result);
                    float_matrix_destroy(&error_matrix);
                }
                i+=10;
                if(i==1031){
                    i=MAX_MAC_SIZE;
                }
                write_to_csv_error(ERROR_FILE, error_array, A_rows,comm_dim,B_cols);
                write_to_csv_latency(LATENCY_FILE, latency,A_rows,comm_dim,B_cols);
                //escribir en archivo
            }
        p++;
        }
    }
    close(memory_file);
    //system("pause");
    return 0;
}