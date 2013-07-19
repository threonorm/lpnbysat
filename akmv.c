#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "limits.h"
#include "omp.h"
#include "mpi.h"

#define CSTVAUDENAY 10000
#define SIZE_32 1 
#define NB_IT 32000*4*SIZE_32

typedef unsigned int ui;

void* print_binary1(unsigned int a,unsigned int k)
{
 if (a==0)
 {int j;for(j=1;j<=k;j++){printf("0");}}
 else{print_binary1(a/2,k-1);printf("%d",a%2);fflush(stdout);} 
}

void* print_binary(unsigned int a,unsigned int k)
{
 if (a==0)
 {int j;for(j=1;j<=k;j++){printf("0\n");}}
 else{print_binary(a/2,k-1);printf("%d\n",a%2);fflush(stdout);} 
}


void* print_mat (unsigned int **mat,int size,int size2)
{
  int i;
  int j;
  for(j=0;j<size;j++)
  { 
    for(i=0;i<size2;i++)
    {
      print_binary1(mat[j][i],32); 
    }
    printf("\n");
  }
  printf("\n");
  fflush(stdout);
}

void* print_vect(unsigned int *a, int size)
{
  int i;
  
  for(i=0;i<size;i++)
  {
    print_binary1(a[i],32); 
  }
  printf("\n");
  fflush(stdout);
}


void* generate_A(ui **mat,int size,int size2)
{
  int i;
  int j;
  for(i=0;i<size;i++)
  {
    for(j=0;j<size2;j++)
    {
      if(rand()%2==0){mat[i][j] = rand();}
      else{mat[i][j]=rand()*2;}
    }
  }
}

void* zero_vect(ui* vect,int size)
{
    int i;
    for(i=0;i<size;i++)
    {
      vect[i] = 0;
    }
}

void* add_noise(ui *v1,ui* v2,float eps,int size)
{
    ui a;
    float x;
    ui rest;
    ui bloc;
    for(a=0;a<size*32;a++)
    {
       rest=a%32;
       bloc=a/32;
       x = (float)rand()/(float)(RAND_MAX);
       if (x<=eps)
       {v2[bloc]=v2[bloc]^(v2[bloc]&(1<<(31-rest)))^(!(v1[bloc]&(1<<(31-rest)))*(1<<(31-rest)));}
       else{v2[bloc]=v2[bloc]^(v2[bloc]&(1<<(31-rest)))^(v1[bloc]&(1<<(31-rest)));}

    }
} 

void* times (ui **mat, ui* vect,ui* result,int size)
{
    int i;
    int j;
    ui bloc,rest;
    ui sca;
    zero_vect(result,SIZE_32);
    for(i=0;i<size*32;i++)
    {
        bloc=i/32;
        rest=i%32;
        sca=0;
        for(j=0;j<SIZE_32;j++)
        {
          int val;
          val = vect[j]&mat[i][j];
          val  = val ^(val>>16);
          val  = val ^(val>>8);
          val  = val ^(val>>4);
          val  = val ^(val>>2);
          val  = val ^(val>>1);
          val  = val &1; 
          sca += val;
        }
        result[bloc]=result[bloc]+(sca%2)*(1<<(31-rest));
    }
}


ui binom(ui n, ui p)
{
  ui binom;
  ui k;			
  if (p > n-p) 
    p = n-p;
  for (binom = 1, k = 1; k <= p; k++)
    binom = (binom * (n-p+k)) / k;
  return binom;
}

void* hamming_bij(ui n, ui w, unsigned long i,ui* temp)
{
  unsigned long res = binom(n,w);
  if(i>=res)
  {
      temp[n/32] = temp[n/32]^(1<<(31-n%32));
      if(n==0){return NULL;}
      hamming_bij(n-1,w-1,i%(binom(n+1,w)-res),temp);
  }
  else
  {
      if(n==0){return NULL;}
      hamming_bij(n-1,w,i,temp);
  }
}

void* generate_hamming(ui** work,ui ndepart,ui n,ui w,ui begin,ui end)
/* Generate all n-words of hamming weight w in work.
 * begin and end are initialized with 0 and binomial(w,n) respectively
 * ndepart=n initially. Actually ndepart is a constant
 */  
{
    ui pos;
    if(w==0||n==0){return NULL;} 
    if(n==w)
    {
       ui j,i;
       for(j=(ndepart-n);j<ndepart;j++)
       {
        for(i=begin;i<end;i++)    
        {
           work[i][j/32]=work[i][j/32]^(1<<(31-j%32));
        }
       } 
       return NULL;
    }
    else
    {
       pos=binom(n,w)-binom(n-1,w);
        if(pos==0){printf("test\n");}
        ui i;
        for(i=begin;i<begin+pos;i++)    
        {
            work[i][(ndepart-n)/32]=work[i][(ndepart-n)/32]^(1<<(31-(ndepart-n)%32));
        }
        generate_hamming(work,ndepart,n-1,w-1,begin,begin+pos);
        generate_hamming(work,ndepart,n-1,w,begin+pos,end);
    }
}


void* xor_btb(unsigned int *a,unsigned int *b,unsigned int *res,int size)
//ICI ON PEUT GAGNER EN SSE
{
  int i;
  for(i=0;i<size;i++)
  {
      res[i]= a[i]^b[i]; 
  }
}

float AKMV_estimator(ui** bici,ui* faibi,ui* faici,ui* c)
/* This function computes an estimation of f^(c), the Walsh-Hadamard
 * transform in c.
 */
{
    float result=0;
    int times=0;
    int j,n;
    for(n=0;n<CSTVAUDENAY;n++)
    {
        times=0;
        for(j=0;j<SIZE_32;j++)
        {
          int val;
          val = c[j] & bici[n][j];
          val  = val ^ (val>>16);
          val  = val ^ (val>>8);
          val  = val ^ (val>>4);
          val  = val ^ (val>>2);
          val  = val ^ (val>>1);
          val  = val &1; 
          times += val;
        }
        if(times%2) {
                  result+= (2*!!(faibi[n/32]&(1<<(31-n%32)))-1)
                  *(2*!!(faici[n/32]&(1<<(31-n%32)))-1)
                  *(-1);
        }
        else{
                 result+= (2*!!(faibi[n/32]&(1<<(31-n%32)))-1)
                 *(2*!!(faici[n/32]&(1<<(31-n%32)))-1);
        }
    }
    result=result/CSTVAUDENAY;
    return result;
}

void* generate_4_AKMV(ui** work,ui* b,ui** bici,ui* faibi,ui* faici)
/*This function makes the datas for akmv estimator*/
{
    int n;
    for(n=0;n<CSTVAUDENAY;n++)
    {
        faici[n/32]=faici[n/32]^((1<<(31- n%32 ))*!!(b[(2*n)/32]&(1<<(31-(2*n)%32))));
        faibi[n/32]=faibi[n/32]^((1<<(31- n%32 ))*!!(b[(2*n+1)/32]&(1<<(31-(2*n+1)%32))));
        xor_btb(work[2*n+1],work[2*n],bici[n],SIZE_32);
    }
}

void* bruteforce_attack(ui* b,ui** samples,ui sizebrute, ui n, ui p)
{
    int i;
    float value;
    //----Malloc time
    
    ui* faibi = malloc(sizeof(ui)*CSTVAUDENAY);
    ui* faici = malloc(sizeof(ui)*CSTVAUDENAY);
    ui* bici1 = malloc(sizeof(ui)*CSTVAUDENAY*SIZE_32);
    long* bici=malloc(sizeof(long)*CSTVAUDENAY);
    ui* temp = malloc(sizeof(ui)*SIZE_32); 
    for(i = 0 ; i < CSTVAUDENAY ; i++)
      bici[i]=(long)&(bici1[SIZE_32*i]);  

    //-----Body
    generate_4_AKMV(samples, b,(ui**)bici,faibi,faici);


    float local;
    local=0;
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
    int compteur=0; 
    //#pragma omp parallel for private(value) 
    for(i = world_rank*sizebrute/world_size; i < (world_rank+1)*sizebrute/world_size; i++)
    {
      zero_vect(temp,SIZE_32);
      hamming_bij(n,p,i,temp);
      value = AKMV_estimator((ui**)bici,faibi,faici,temp);
      if(value>=0.01){compteur++;
			}
    }
    printf("%d",compteur);
    //-----Free time
    free(faici);
    free(faibi);
    free(bici); 
    free(bici1);
}

int main()
{
  int i;
  int nberreur = 3;
  omp_set_num_threads(32);
  srand(time(NULL));
  
  unsigned long t = binom(32*SIZE_32,nberreur);
  //---Malloc time
  
  long* work = malloc(sizeof(long)*NB_IT);
  ui* work1= malloc(sizeof(ui)*NB_IT*(SIZE_32));
  long* hamming = malloc(sizeof(long)*t);
  ui* hamming1 = malloc(sizeof(ui)*SIZE_32*t);
  ui* s = malloc(sizeof(ui)*SIZE_32);
  ui* b = malloc(sizeof(ui)*NB_IT/32); 
  
  for(i=0;i<t;i++)
  {
     hamming[i] = (long)&(hamming1[SIZE_32*i]); 
  }

  for(i=0;i<NB_IT;i++)
  {
     work[i] = (long)&(work1[SIZE_32*i]); 
  }
  //----MPI gestion
  
  MPI_Init(NULL,NULL);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
   
  
  //----Body
  
  if (world_rank == 0)
  {
    //Generate data  
    generate_A((ui**)work,NB_IT,SIZE_32);
    zero_vect(b,NB_IT/32);
    s[0] = (1<<31) + (1<<14) + (1<<1);   
    times((ui**)work, s, b, NB_IT/32);
    add_noise(b, b, 0.45, NB_IT/32);
  }
 
  //Broadcast 
  MPI_Bcast(b, NB_IT / 32, MPI_INT, 0, MPI_COMM_WORLD); 
  int recup;
  for(recup = 0; recup < NB_IT ; recup++)
  {
    MPI_Bcast(&work1[SIZE_32*recup], SIZE_32, MPI_INT, 0, MPI_COMM_WORLD);
  }

  

  //Attack
  bruteforce_attack((ui*)b, (ui**)work, t, SIZE_32*32, nberreur);

  
  //----Free time
  free(work1);
  free(work);
  free(hamming1);
  free(hamming);
  free(s);
  free(b);
  //-------- 
  return 0;
}
