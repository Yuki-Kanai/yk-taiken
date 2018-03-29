/*
   To-do: optimization
   handle multiple meshes
   */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define cut_num 5
#define LP_root {0.148874338981631, 0.433395394129247, .679409568299024, 0.865063366688985, 0.973906528517172}
#define LP_weight {0.29552422471475,  0.269266719309996, 0.21908636251598, 0.149451349150581,  0.066671344308688}

//x0 of element
double get_x(int elem_num){
    return (elem_num % cut_num) / (double)cut_num;
}

double get_y(int elem_num){
    return (elem_num / cut_num) / (double) cut_num;
}

//fill array with 0
void reset_matrix_double(double* mat, int num){
    int i;
    for(i = 0; i!=num; i++)
        mat[i] = 0;
}

void mult_matrix(const double* A, const double* B, double* C, int ma, int mb, int nb){
    int i,j,k;
    reset_matrix_double(C, ma*nb);
    for(i=0; i!= ma; i++){
        for(j=0; j!= nb; j++){
            for(k=0; k!= mb; k++){
                C[i*nb+j] += A[i*mb+k] * B[k*nb+j];
            }}}
}

//Product of two shape functions N1(x,y)*N1, N1*N2, ..., Return det J
double SF_generator(double x, double y, int elem_num, double* out_array, int is_NN){
    double g1 = (x- get_x(elem_num))*2 * cut_num - 1;
    double g2 = (y- get_y(elem_num))*2 * cut_num - 1;
    //dN1/dg1, dN1/dg2, dN2/dg1, ...
    double dNg_array[8] = {-(1-g2)/4, -(1-g1)/4, (1-g2)/4, -(1+g1)/4, (1+g2)/4, (1+g1)/4, -(1+g2)/4, (1-g1)/4};
    int i, j;
    double J11 = 2* cut_num, J12 = 0,
           J21 = 0, J22 = 2* cut_num;
    double detJ = J11*J22 - J21*J12;
    //To-do stderr if Ji==0
    if(is_NN){
        double temp[4];
        temp[0] = (1-g1)*(1-g2)/4;
        temp[1] = (1+g1)*(1-g2)/4;
        temp[2] = (1+g1)*(1+g2)/4;
        temp[3] = (1-g1)*(1+g2)/4;
        for(i=0; i!=4; i++)
            for(j=0; j!=4; j++)
                out_array[i*4+j] = temp[i]*temp[j];
    }else{
        for(i=0; i!=4; i++)
            for(j=0; j!=4; j++){
                out_array[i*4+j] = (J11 * dNg_array[i*2] + J12 * dNg_array[i*2+1]) * (J11 * dNg_array[j*2] + J12 * dNg_array[j*2+1])+
                    (J21 * dNg_array[i*2] + J22 * dNg_array[i*2+1]) * (J21 * dNg_array[j*2] + J22 * dNg_array[j*2+1]);
            }
    }
    return detJ;
}


//output 16 Integrations throughout 2D Isoparametric square of products of two shape functions
void get_localNN(double* localNN, int elem, int is_NN){
    int i,j,k,l,m;
    double x0 = get_x(elem), y0 = get_y(elem);
    double LP_roots[5] = LP_root, LP_weights[5] = LP_weight;
    double Jacobs[4];
    double temp[4][16];//for each 2*2 GL integration field
    reset_matrix_double(localNN, 16);
    for(i=0; i!=5; i++){
        for(j=0; j!=5; j++){
            Jacobs[0] = SF_generator(x0+(1+LP_roots[i])/(double)2/ cut_num, y0+(1+LP_roots[j])/(double)2/ cut_num, elem, temp[0], is_NN);
            Jacobs[1] = SF_generator(x0+(1-LP_roots[i])/(double)2/ cut_num, y0+(1+LP_roots[j])/(double)2/ cut_num, elem, temp[1], is_NN);
            Jacobs[2] = SF_generator(x0+(1-LP_roots[i])/(double)2/ cut_num, y0+(1-LP_roots[j])/(double)2/ cut_num, elem, temp[2], is_NN);
            Jacobs[3] = SF_generator(x0+(1+LP_roots[i])/(double)2/ cut_num, y0+(1-LP_roots[j])/(double)2/ cut_num, elem, temp[3], is_NN);
            for(k=0; k!=4; k++){
                for(l=0; l!=4; l++){
                    for(m=0; m!=4; m++)
                        localNN[k*4+l]+=temp[m][k*4+l]*LP_weights[i]*LP_weights[j] * Jacobs[m];
                }}}}
}

int get_global_index(int elem, int node){
    switch(node){
        case 0:
            return (elem/cut_num) * (cut_num+1) +  elem%cut_num;
        case 1:
            return (elem/cut_num) * (cut_num+1) +  elem%cut_num + 1;
        case 2:
            return (elem/cut_num + 1) * (cut_num+1) +  elem%cut_num + 1;
        case 3:
            return (elem/cut_num + 1) * (cut_num+1) +  elem%cut_num;
    }
    return -1;
}

void flip_row(double* A, int i, int j, int size){
    double* array;
    array = (double *)malloc(sizeof(double)*size);
    int k=0;
    for(k=0; k!=size; k++)
        array[k] = A[i*size+k];
    for(k=0; k!=size; k++)
        A[i*size+k] = A[j*size+k];
    for(k=0; k!=size; k++)
        A[j*size+k] = array[k];
    free(array);
}

//LU decomposition of size*size square matrix,C: input -> upper triangle ,L: lower triangle 、P: Permutation matrix for Pivot Selection、S:scaling
void LU_factorize(double* C, double* L, double* P, double* S, int size){
    printf("Start LU factorization\n");
    reset_matrix_double(L, size*size);
    reset_matrix_double(P, size*size);
    reset_matrix_double(S, size*size);
    int i,j,k,l, c_max_i;
    double c_max, temp;
    //scale
    for(i=0; i!=size; i++){
        c_max = fabs(C[i*size]);
        c_max_i = 0;
        for(j=1; j!=size; j++){
            if(c_max < fabs(C[i*size+j])){
                c_max = fabs(C[i*size+j]);
                c_max_i = j;
            }
        }
        S[i*size+i] = 1/c_max;
        for(j=0; j!=size; j++)
            C[i*size+j] /=c_max;
    }
    printf("End Scaling\n");
    for(i=0; i!=size; i++){
        P[i*size+i] = 1;
        L[i*size+i] = 1;
    }
    l=0;
    printf("Start Gaussian Elimination\n");
    for(i=0; i!= size-1; i++){
        //Pivot selection
        c_max = fabs(C[i*size+i]);
        c_max_i = i;
        for(j=i+1; j!=size; j++){
            if(fabs(C[j*size + i]) > c_max){
                c_max =fabs(C[j*size+i]);
                c_max_i = j;
            }
        }
        flip_row(P, i, c_max_i, size);//stores past exchanges
        flip_row(C, i, c_max_i, size);//A(i) before gaussian delete
        flip_row(L, i, c_max_i, size);//

        for(j=i+1; j!=size; j++){
            temp = C[j*size+i]/C[i*size+i];
            L[j*size+i] = temp;
            C[j*size+i] = 0;
            for(k=i+1; k!=size; k++)
                C[j*size+k] -= -temp * C[i*size+k];
        }

        if((i*10)/size > l){
            printf("\rGauss Elimination on progress:%d0%% completed ",l);
            fflush(stdout);
            l++;
        }
    }
    printf("End LU factorization\n\n");
}

//forward substition of LU decomposition :Lc=b
void forward_elimination(const double* L, double* c, const double * b, int n){
    int i,j;
    double temp;
    for(i=0; i!=n; i++){
        temp = b[i];
        for(j=0; j!=i; j++)
            temp -= b[j] * L[i*n+j];
        c[i] = temp/L[i*n+i];
    }
}

//backward substition of LU decomposition Ux=c
void backward_substitution(const double* U, double* x, const double* c, int n){
    int i,j;
    double temp;
    for(i=n-1; i!=-1; i--){
        temp = c[i];
        for(j=n-1; j!= i; j--)
            temp -= c[j]*U[i*n+j];
        x[i] = temp/U[i*n+i];
    }
}

void print_status(FILE *fp, const double* phi, int size, int i, double h){
    int j;
    for(j=0; j!=size; j++){
        fprintf(fp, "%d %f %f\n", j, i*h, phi[j]);
    }
}

//insert array must have four times the maximum number of nodes of each type
//in:node inside the submesh, bound, node: inside but shares lattices with other submeshes, out: outside but shares lattices with this submesh 
void set_in_global(int *in_global_i, int* in_nums){
    int i,j;
    int cut_size;
    if(cut_num%2==1){
        cut_size = (cut_num+1)/2;
        for(i=0; i!=cut_size; i++){
            for(j=0; j!=cut_size; j++){
                in_global_i[i*cut_size+j] = i*(cut_num+1)+j;
                in_global_i[(cut_size*cut_size)+i*cut_size+j] = i*(cut_num+1)+j+cut_size;
                in_global_i[(cut_size*cut_size)*2+i*cut_size+j] = (i+cut_size)*(cut_num+1)+j;
                in_global_i[3*(cut_size*cut_size)+i*cut_size+j] = (i+cut_size)*cut_num+j+cut_size;
            }}
        for(i=0; i!=4; i++)
            in_nums[i]=cut_size*cut_size;
    }
    else{
            cut_size = cut_num/2;//value would be small end of rectangle
            for(i=0; i!=cut_size; i++){//short
                for(j=0; j!=cut_size+1; j++){//long
                    in_global_i[i*(cut_size+1)+j] = i*(cut_num+1)+j;
                    in_global_i[i*(cut_size+1)+j+(cut_size+1)*cut_size] = i+cut_size+1 + (cut_num+1)*j;
                    in_global_i[i*(cut_size+1)+j+(cut_size+1)*cut_size*2] = (i+cut_size+1)*(cut_num+1) + cut_size+j;
                    in_global_i[i*(cut_size+1)+j+(cut_size+1)*cut_size*3] = i + (j+cut_size)*(cut_num+1);
                }}
        in_global_i[(cut_size+1)*cut_size*4] = cut_size + cut_size*(cut_num+1);

        for(i=0;i!=3;i++){
            in_nums[i] = cut_size*(cut_size+1);
        }
        in_nums[3] = cut_size*(cut_size+1)+1;
    }
}

void divide_mesh(int* in_global_i, int *bound_global_i, int *out_global_i, int* in_nums, int* bound_nums, int* out_nums,
        int* elems, int*elem_nums){
    int i,j;
    int cut_size;
 
    set_in_global(in_global_i, in_nums);
    set_bound_global(bound_global_i, bound_nums);
    if(cut_num%2==1){
        cut_size = (cut_num+1)/2;
        for(i=0; i!=cut_size; i++){
            bound_global_i[i] = (cut_size-1)*(cut_num+1)+i;
            bound_global_i[i+(cut_size*2-1)] = (cut_size-1)*(cut_num+1)+i+cut_size;
            bound_global_i[i+(cut_size*2-1)*2] = (cut_size)*(cut_num+1)+i;
            bound_global_i[i+(cut_size*2-1)*3] = (cut_size)*(cut_num+1)+i+cut_size;
            if(i!=cut_size-1){
                bound_global_i[i+cut_size] = cut_size-1+(cut_num+1)*i;
                bound_global_i[i+cut_size+(cut_size*2-1)*2] = cut_size + (cut_num+1)*i;
                bound_global_i[i+cut_size+(cut_size*2-1)*3] = cut_size + (cut_num+1)*(i+cut_num);
                bound_global_i[i+cut_size+(cut_size*2-1)*4] = cut_size -1 + (cut_num+1)*(i+cut_num);
            }
            out_global_i[i] = cut_size*(cut_num+1) + i;
            out_global_i[i+cut_size] = cut_size + (cut_num+1)*i;
            out_global_i[i+cut_size*2+1] = cut_size*(cut_num+1) + cut_size + i;
            out_global_i[i+cut_size*2+1+cut_size] = cut_size - 1 + (cut_num+1)*i;
            out_global_i[i+(cut_size*2+1)*2] = (cut_size-1)*(cut_num+1)+i;
            out_global_i[i+(cut_size*2+1)*2+cut_size] = cut_size - 1 + (cut_num+1)*(i+cut_size-1);
            out_global_i[i+(cut_size*2+1)*3] = (cut_size-1)*(cut_num+1)+i+cut_size;
            out_global_i[i+(cut_size*2+1)*4] = cut_size + (cut_num+1)*(i+cut_size-1);
            out_global_i[cut_size*2]= cut_size + cut_size*(cut_num+1);
        }
        out_global_i[(cut_size*2+1)*2-1] = cut_size-1 + cut_size*(cut_num+1);
        out_global_i[(cut_size*2+1)*3-1] = cut_size-1+(cut_size-1)*(cut_num+1);
        out_global_i[(cut_size*2+1)*4-1] = cut_size + (cut_size-1)*(cut_num+1);
        for(i=0; i!=4; i++){
            bound_nums[i]=cut_size*2-1;
            out_nums[i]= cut_size*2+1;
        }
        for(i=0; i!=(cut_size+1);i++)
            for(j=0; j!=(cut_size+1); j++){
                elems[i*cut_size+j] = i*cut_num+j;
                elems[i*cut_size+j+(cut_size+1)*(cut_size+1)] = i*cut_num+j+cut_size;
                elems[i*cut_size+j+(cut_size+1)*(cut_size+1)*2] = (cut_num-1-i)*(cut_num)+j;
                elems[i*cut_size+j+(cut_size+1)*(cut_size+1)*3] = (cut_num-i)*(cut_num)-j-1;
            }
        for(i=0; i!=4; i++)
            elem_nums[i]=cut_size*cut_size;
    }
    else{
        cut_size = cut_num/2;//value would be small end of rectangle
        for(i=0; i!=cut_size; i++){//short
            for(j=0; j!=cut_size+1; j++){//long
                if(i==0){
                    bound_global_i[j] = (cut_num+1)*(cut_size-1)+j;
                    bound_global_i[j+cut_size*2] = (cut_num+1)*j + cut_size+1;
                    bound_global_i[j+cut_size*2*2] = (cut_num+1)*(cut_size+1) + cut_size+j;
                    bound_global_i[j+cut_size*2*3] = (cut_num+1)*cut_size + j;

                    out_global_i[i] = (cut_size)*(cut_num+1)+i;
                    out_global_i[i+(cut_size*2+2)] = cut_size + i+(cut_num+1);
                    out_global_i[i+(cut_size*2+2)*2] = (cut_size)*(cut_num+1) + (cut_num-i);
                    out_global_i[i+(cut_size*2+2)*3] = (cut_size-1)*(cut_num+1) + i;
                }
            }
            bound_global_i[cut_size+i] = cut_size+i*(cut_num+1);
            bound_global_i[i+cut_size*3] = (cut_num+1)*cut_size + cut_size+1+i;
            bound_global_i[i+cut_size*5] = (cut_num+1)*(cut_size+1+i) + cut_size;
            bound_global_i[i+cut_size*7] = (cut_num+1)*cut_size + i;

            out_global_i[i+cut_size+2] = cut_size+1 + (cut_num+1)*i;
            out_global_i[(cut_size*2+2)*2-1-i] = (cut_num+1)*cut_size - 1 -i;
            out_global_i[(cut_size*2+2)*3-1-i] = cut_size-1 + (cut_num+1)*(cut_num-i);
            out_global_i[(cut_size*2+2)*4-1-i] = cut_size + (cut_num+1)*(cut_num-i);
        }
        bound_global_i[cut_size*8] = cut_size*(cut_num+1+1);
        for(i=0; i!=3; i++){
            out_global_i[(cut_size*2+2)*4+i] = (cut_num+1)*(cut_num-i) + cut_size;
        }
        for(i=0;i!=3;i++){
            bound_nums[i] = cut_size*2;
            out_nums[i] = (cut_size+1)*2;
        }
        bound_nums[3] = cut_size*2+1;
        out_nums[3] = (cut_size+1)*2 +2;
        for(i=0; i!=cut_size; i++)
            for(j=0; j!=(cut_size+1); j++){
                elems[i*(cut_size+1)+j] = j + i*cut_num;
                elems[(cut_size+i)*(cut_size+1)+j] = cut_num-1-i+j*cut_num;
                elems[(cut_size*2+i)*(cut_size+1)+j] = cut_num*(cut_num-i)-j-1;
                elems[(cut_size*3+i)*(cut_size+1)+j] = cut_num*(cut_num-j-1)+ i;
            }
        elems[cut_size*(cut_size+1)*4]=cut_num*(cut_num-cut_size)+ cut_size;
        elems[cut_size*(cut_size+1)*4+1]=cut_num*(cut_num-cut_size-1)+ cut_size;
        for(i=0; i!=3;i++)
            elem_nums[i]=cut_size*(cut_size+1);
        elem_nums[3]=cut_size*(cut_size+1)+2;
    } 
}

void set_bound(int* fix_bound_i, int* fix_bound_num){
    *fix_bound_num = (cut_num+1)*2;
    for(i=0;i!=cut_num+1;i++){
        fix_bound[i] = 100;
        fix_bound_i[i] = i*(cut_num+1);
    }
    for(i=0; i!= (cut_num+1);i++){
        fix_bound[i+cut_num+1] = 0;
        fix_bound_i[i+cut_num+1] = (i+1)*(cut_num+1) -1;
    }
}

int get_local_index(int global_i, int* local_nodes, int local_node_tot){
    int local_i = -1;
    int i;
    for(i=0; i<= local_node_tot; i++){
        if(local_nodes[i]=global_i){
            local_i=i;
            break;
        }}
    return local_i
}

void set_KC(double* K, double* C, int* subfield_i, int sub_field_tot, int* node_i, int node_tot,
        int* fix_bound_i, int fix_bound_num){
    double temp_localNN[16], temp_localdNN[16];

    reset_matrix_double(NN, node_tot*node_tot);
    reset_matrix_double(dNN, node_tot*node_tot);

    k=0;
    printf("Calculate element stiffness equation for each elements\n");
    for(i=0; i!=sub_field_tot; i++){
        get_localNN(temp_localNN, i, 1);//16 products
        get_localNN(temp_localdNN, i, 0);
        for(j=0; j!=4; j++){
            for(k=0; k!=4; k++){
                NN[get_global_index(i,j)* size+ get_global_index(i,k)] += temp_localNN[j*4+k];
                dNN[get_global_index(i,j)* size+ get_global_index(i,k)] += temp_localdNN[j*4+k]*D; //warning: multiply D
            }
        }
        if((i*10)/(sub_field_tot>k)){
            k++;
            printf("\rCalculating ES equation：%d0%%",k);
            fflush(stdout);
        }
    }
    
    int index;
    for(i=0; i!= fix_bound_num; i++){
        for(j=0; j!=size; j++){
            index = get_local_index(fix_bound_i[i], node_i, node_tot);
            if(index >= 0)
                NN[index] = 0;
        }
    }
    for(i=0; i!= fix_bound_num; i++){
        index = get_local_index(fix_bound_i[i], node_i, node_tot);
        if(intdex >= 0)
            NN[index] = 1;
    }
    printf("\nFinished calculating ES Equation\n");
}

void print2stdout(double* phi, int size, int i){
    printf("\n\nfinished calculating until i=%d\n\n",i);
    printf("phi[i]=");
    for(j=0; j!=size; j++)
        printf("%.2e ",phi[size]);
    printf("\n\n");

}

int main(int argc, char* argv[]){
    FILE* out;
    int i,j,k;
    double h = 1e-6;//time interval
    int itmax = 1e6;
    double D = 1;
    double* NN, *dNN, *C, *L, *P, *S, *temps;
    int *in_global_i, *bound_global_i, *out_global_i, in_nums[4], bound_nums[4], out_nums[4], *elems, elem_nums[4];

    int in_num, bound_num, out_num;
    if(cut_num%2==1){
        in_num=(cut_num+1)*(cut_num+1); bound_num=((cut_num/2)*2+1)*4; out_num=((cut_num/2)*2+3)*4;
    }else{
        in_num=(cut_num+1)*(cut_num+1); bound_num=cut_num*4+1; out_num=(cut_num+2)*4+2;
    }
    in_global_i=(int*)malloc(sizeof(int)*in_num); bound_global_i=(int*)malloc(sizeof(int)*bound_num);
    out_global_i=(int*)malloc(sizeof(int)*out_num); elems=(int*)malloc(sizeof(int)*in_num);
    divide_mesh(in_global_i, bound_global_i, out_global_i, in_nums, bound_nums, out_nums, elems, elem_nums);


    double* fix_bound;
    int fix_bound_num; 
    fix_bound = (double*)malloc(sizeof(double)*fix_bound_num);
    int* fix_bound_i;
    fix_bound_i=(int*)malloc(sizeof(int)*fix_bound_num);
    set_bound(fix_bound, &fix_bound_num);
    
    int rank, proc;
    MPI_INIT(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    
    int size = in_num[rank]+out_num[rank];
    NN = (double*)malloc(sizeof(double)*size*size);
    dNN = (double*)malloc(sizeof(double)*size*size);
    C = (double*)malloc(sizeof(double)*size*size);
    L = (double*)malloc(sizeof(double)*size*size);
    P = (double*)malloc(sizeof(double)*size*size);
    S = (double*)malloc(sizeof(double)*size*size);
    temps = (double*)malloc(sizeof(double)*size*size);
    double temp_localNN[16], temp_localdNN[16];

    void set_KC(K, C, subfield_i, sub_field_tot, node_i, node_tot, fix_bound_i, fix_bound_num);
    LU_factorize(NN, L, P, S, size);

    double* b, *x, *phi, *temp;
    b = (double*)malloc(sizeof(double)*size);
    x = (double*)malloc(sizeof(double)*size);
    phi = (double*)malloc(sizeof(double)*size);
    temp = (double*)malloc(sizeof(double)*size);
    if(
    out = fopen("result.data","w");


    reset_matrix_double(phi, size);//set to 0degrees for initial condition
    print_status(out, phi, size, 0, h);

    for(i=0; i!= itmax; i++){
        mult_matrix(dNN, phi, b, size, size, 1);
        for(j=0; j!=size; j++)
            b[j] = -b[j];
        //to set dphi/dt to 0 for Dirichlet boundary
        for(j=0; j!=fix_bound_num; j++)
            b[fix_bound_i[j]] = 0;

        mult_matrix(S, b, temp, size, size, 1);
        mult_matrix(P, temp, b, size, size, 1);

        forward_elimination(L, temp, b, size);
        backward_substitution(NN, x, temp, size);

        for(j=0; j!=size; j++)
            phi[j] += x[j]*h;
        //For this case, the boundary is set by dphi/dt, so no change is required for vector b
        for(j=0; j!=fix_bound_num; j++){
            phi[fix_bound_i[j]] = fix_bound[j];
        }
        print_status(out, phi, size, i, h);
        if(i%(itmax/10)==0) 
            print2stdout(phi, size);
    }
    MPI_Finalize();
    printf("\n\n!!!Finished!!!\n\n");
    free(in_global_i); free(bound_global_i); free(out_global_i); free(elems);
    free(NN); free(dNN); free(C); free(L); free(P); free(S); free(temps);
    free(b); free(x); free(phi); free(temp);
    free(fix_bound); free(fix_bound_i);
    fclose(out);
    return 0;
}

