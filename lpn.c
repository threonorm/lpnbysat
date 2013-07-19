#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "limits.h"
#include "omp.h"
#define CSTVAUDENAY 10000
#define SIZE_32 1
#define NB_IT 70000*SIZE_32

typedef unsigned int ui;

//Gauss + fcts bases lpn etc...


ui max(ui *tab,int size)
{
  int i;
  ui loc=0;
  ui res;
  for(i=0;i<32*size;i++)
  {
    if(tab[i]>=loc)
      {
          loc=tab[i];
          res=i;
      } 
  }
  return res;
}


void* zero_vect(unsigned int* vect,int size)
{
    int i;
    for(i=0;i<size;i++)
    {
      vect[i] = 0;
    }
}


ui pow1 (ui a, ui b)
{
  if(b==0){return 1;}
  if(b==1){return a;}
  else{if(b%2==0){return pow1(a*a,b/2);}
        else{return (a*pow1(a*a,b/2));}}  
}


void* zeroA(unsigned long **mat, int size,int size2)
{
  int i;
  int j;
  for(i=0;i<size;i++)
  {
    for(j=0;j<size2;j++)
    {
      mat[i][j] = 0;
    }
  }
}


void* zero_A(unsigned int **mat, int size,int size2)
{
  int i;
  int j;
  for(i=0;i<size;i++)
  {
    for(j=0;j<size2;j++)
    {
      mat[i][j] = 0;
    }
  }
}

void* generate_A(unsigned int **mat,int size,int size2)
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
 
void* transpose_A(unsigned int **mat,int size,int size2)
{
  int i;
  int j;
  for(i=0;i<32*size;i++)
  {
    for(j=i+1;j<32*size2;j++)
    {
      ui tempij;
      tempij=!!(mat[i][j/32]&(1<<(31-j%32)));
      mat[i][j/32]=mat[i][j/32]^(mat[i][j/32]&(1<<(31-j%32)))^(!!(mat[j][i/32]&(1<<(31-i%32)))*(1<<(31-j%32)));
      mat[j][i/32]=mat[j][i/32]^(mat[j][i/32]&(1<<(31-i%32)))^(tempij*(1<<(31-i%32)));
    }
  }
}

void* generate_vect(unsigned int* vect,int size)
{
    int i;
    for(i=0;i<size;i++)
    {
      if(rand()%2==0){vect[i] = rand();}
      else{vect[i]=rand()*2;}
    }
}

void* xor_btb(unsigned int *a,unsigned int *b,unsigned int *res,int size)
{
  int i;
  for(i=0;i<size;i++)
  {
      res[i]= a[i]^b[i]; 
  }
}

void* print_binary(unsigned int a,unsigned int k)
{
 if (a==0)
 {int j;for(j=1;j<=k;j++){printf("0");}}
 else{ print_binary(a/2,k-1);printf("%d",a%2);fflush(stdout);} 
}

typedef
struct cons
{
  ui head;
  struct cons *tail;  
}
cons;

typedef cons* list;



void* exchange(ui **tps0,ui **tps1, ui i , ui j,int size_col)
{
  int l;
  for(l=0;l<size_col;l++)
  {
    int temp;
    temp=tps1[i][l];
    tps1[i][l]= tps0[j][l];
    tps1[j][l]= temp;     
  }
  return NULL; 
}

void* mix(ui **tps0,ui **tps1, ui i, ui j,int size_col)
//matrice carré : 32 fois plus de lignes que de colonnes
{
    xor_btb(tps0[i],tps0[j],tps1[i],size_col);
}


void* print_vect(unsigned int *a, int size)
{
  int i;
  
  for(i=0;i<size;i++)
  {
    print_binary(a[i],32); 
  }
  printf("\n");
  fflush(stdout);
}

void* printmat (unsigned long **mat,int size,int size2)
{
  int i;
  int j;
  for(j=0;j<size;j++)
  { 
    for(i=0;i<size2;i++)
    {
      print_binary(mat[j][i],32); 
    }
    printf("\n");
  }
  printf("\n");
  fflush(stdout);
}

void* print_mat (unsigned int **mat,int size,int size2)
{
  int i;
  int j;
  for(j=0;j<size;j++)
  { 
    for(i=0;i<size2;i++)
    {
      print_binary(mat[j][i],32); 
    }
    printf("\n");
  }
  printf("\n");
  fflush(stdout);
}
 
ui hamming_weight_int(ui a)
{
 ui weight=0;
 for(;a>0; a &=a-1) weight++;
 return weight;
}

ui hamming_weight(ui *a,int size)
{
  ui weight=0;
  int i;
  for(i=0;i<size;i++)
  {
    weight+= hamming_weight_int(a[i]);  
  }
return weight;
}

void* incr(ui vec[SIZE_32])
{
  ui i= SIZE_32-1;
  while(vec[i]==UINT_MAX)
  {i--;}
  vec[i]++;
  if(vec[i]%1000000==0)
  {printf("%u \n",vec[i]); fflush(stdout);}
}

void* times1 (unsigned int **mat, unsigned int* vect,unsigned int* result,int size)
{
    int i;
    int j;
    ui bloc,rest;
    ui sca;
    zero_vect(result,size);
    for(i=0;i<size*32;i++)
    {
        bloc=i/32;
        rest=i%32;
        sca=0;
        for(j=0;j<size;j++)
        {
            //BITTRICK
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
        result[bloc]=result[bloc]+(sca%2)*(1<<(31-rest)/*pow1(2,31-rest)*/);
    }
}


void* times (unsigned int **mat, unsigned int* vect,unsigned int* result,int size)
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
            //BITTRICK
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
        result[bloc]=result[bloc]+(sca%2)*(1<<(31-rest)/*pow1(2,31-rest)*/);
    }
}


ui* exhaustive(ui **mat,ui vec[NB_IT],ui result[SIZE_32])
{
    ui min_bad_values;
    ui s[SIZE_32];
    ui work1[NB_IT];
    ui work2[NB_IT];
    ui ham_work;
    ui ham_resul=UINT_MAX;
    zero_vect(s,SIZE_32);
    while(s[0] != UINT_MAX)
    {
      times(mat,s,work1,SIZE_32);
      xor_btb(work1, vec, work2,NB_IT);
      ham_work = hamming_weight(work2,NB_IT);
      if(ham_work< ham_resul)
        {*result=*s;ham_resul=ham_work;}   
      incr(s);
    }
    return result;   
}

void* randomize(ui **matrix, ui *b, ui **hist0, ui size)
{
  ui i;
  ui sw;
  ui *tmp;
  ui blocc1,restc1;
  for(i=0; i<size;i++)
  {
    sw=rand()%(size);
    //tmp=matrix[i];
    xor_btb(matrix[i],matrix[sw],matrix[i],SIZE_32);
        //matrix[i]=matrix[sw];
    //matrix[sw]=tmp;
    mix(hist0,hist0,i,sw,size/32); 
    ui *tmp;
    ui tmp1;
    blocc1=(i)/32;
    restc1=(i)%32;
    tmp1=!!(b[blocc1]&(1<<(31-restc1))); 
    b[blocc1]=b[blocc1]^((!!(b[sw/32]&(1<<(31-sw%32))))*(1<<(31-restc1)));
  /*  b[sw/32]=(tmp1*(1<<(31-sw%32)))^b[sw/32]^(b[sw/32]&((1<<(31-sw%32))));*/
  }
}


void* randomize_1(ui **matrix, ui *b, ui size)
{
  ui i;
  ui sw;
  ui *tmp;
  ui blocc1,restc1;
  for(i=0; i<32*size;i++)
  {
    sw=rand()%(32*size);
    tmp=matrix[i];
    matrix[i]=matrix[sw];
    matrix[sw]=tmp;
    ui *tmp;
    ui tmp1;
    blocc1=(i)/32;
    restc1=(i)%32;
    tmp1=!!(b[blocc1]&(1<<(31-restc1))); 
    b[blocc1]=b[blocc1]^(b[blocc1]&(1<<(31-restc1)))^((!!(b[sw/32]&(1<<(31-sw%32))))*(1<<(31-restc1)));
    b[sw/32]=(tmp1*(1<<(31-sw%32)))^b[sw/32]^(b[sw/32]&((1<<(31-sw%32))));
  }
}


int gauss_test(ui **matrix, ui *b,ui **hist0,ui size)
{
//Gaussian elimination
 ui i,j;
 ui blocc, restc;
 ui k=0;
 ui blocc2,restc2,blocc1,restc1;
 ui decal=0;
 for(i=0;i<SIZE_32*32;i++)
 { 
       
        if(k==32*size){decal++;}
        blocc=(i)/32;
        restc=(i)%32;
        k=i+1-decal;
        if(!(matrix[i-decal][blocc]&(1<<(31-restc)/*pow1(2,31-restc)*/)))
        {
            while(k<32*size)
            {
              if((matrix[k][blocc]&(1<<(31-restc)/* pow1(2,31-restc)*/)))
              {
                  exchange(hist0,hist0,i-decal,k,size); 
                  ui *tmp;
                  ui tmp1;
                  tmp=matrix[i-decal];
                  matrix[i-decal]=matrix[k];
                  matrix[k]=tmp;

                  //TODO: Suite à modifier
                  blocc1=(i-decal)/32;
                  restc1=(i-decal)%32;
                  tmp1=!!(b[blocc1]&(1<<(31-restc1)/*pow1(2,31-restc1)*/)); 
                  b[blocc1]=b[blocc1]^(b[blocc1]&(1<<(31-restc1)/*pow1(2,31-restc1)*/))^((!!(b[k/32]&(1<<(31-k%32))/*pow1(2,31-k%32)*/))*(1<<(31-restc1))/*pow1(2,31-restc1)*/);
                  b[k/32]=(tmp1*(1<<(31-k%32))/*pow1(2,31-k%32)*/)^b[k/32]^(b[k/32]&((1<<(31-k%32))/*pow1(2,31- k%32)*/));
                  k=32*size;
              }
              k++;
            }
        }
        if(k!=32*size)
        {
          for(j=0;j<32*size;j++)
          {
             if(j!=i-decal)
             {
                int l;
                if(matrix[j][blocc]&( 1<<(31-restc)/*pow1(2,31-restc)*/))
                {
                  mix(hist0,hist0,j,i-decal,size);
                  for(l=0;l<SIZE_32;l++)
                  {
                    matrix[j][l]= matrix[j][l] ^ matrix[i-decal][l];    
                  }
                  blocc2=(i-decal)/32;
                  restc2=(i-decal)%32;
                  b[j/32]=b[j/32]^(!!(b[blocc2]&(1<<(31-restc2))/*pow1(2,31-restc2)*/)*(1<<(31-j%32)/*pow1(2,31-j%32)*/));
                }  
             }
          }
        }
 }
 ui kl=size*32-1;
 int solu=1;
 int count=0;
 ui *erreur;
 erreur=malloc(sizeof(ui)*size*32);
 int init;
 for(init=0;init<size*32;init++){erreur[init]=0;}
 while(!hamming_weight(matrix[kl],SIZE_32))
 {
     
      if(!!(b[kl/32]&(1<<(31-kl%32))) != 0)
      {
        int p;
        count++;
        for(p=0;p<size*32;p++)
        {
          init=hist0[kl][p/32]&( 1<<(31-p%32)); //SURPRENANT
          if(init)
          {
            erreur[p]++;
          }
        }
        //print_vect(hist0[kl],size);  
      } 
      kl--;
 }

 int final=max(erreur,size);
return kl;
}

int gauss2(ui **matrix, ui *b,ui size,ui larg)
{
//Gaussian elimination
 ui i,j;
 ui blocc, restc;
 ui k=0;
 ui blocc2,restc2,blocc1,restc1;
 ui decal=0;
 for(i=0;i<larg*32;i++)
 {
      print_mat(matrix,32*size,larg); 
        if(k==32*size){decal++;}
        blocc=(i)/32;
        restc=(i)%32;
        k=i+1-decal;
        if(!(matrix[i-decal][blocc]&(1<<(31-restc)/*pow1(2,31-restc)*/)))
        {
            while(k<32*size)
            {
          
              if((matrix[k][blocc]&(1<<(31-restc)/* pow1(2,31-restc)*/)))
              {
                  ui *tmp;
                  ui tmp1;
                  tmp=matrix[i-decal];
                  matrix[i-decal]=matrix[k];
                  matrix[k]=tmp;

                  //TODO: Suite à modifier
                  blocc1=(i-decal)/32;
                  restc1=(i-decal)%32;
                 // tmp1=!!(b[blocc1]&(1<<(31-restc1)/*pow1(2,31-restc1)*/)); 
                 // b[blocc1]=b[blocc1]^(b[blocc1]&(1<<(31-restc1)/*pow1(2,31-restc1)*/))^((!!(b[k/32]&(1<<(31-k%32))/*pow1(2,31-k%32)*/))*(1<<(31-restc1))/*pow1(2,31-restc1)*/);
//                  b[k/32]=(tmp1*(1<<(31-k%32))/*pow1(2,31-k%32)*/)^b[k/32]^(b[k/32]&((1<<(31-k%32))/*pow1(2,31- k%32)*/));
                  k=32*size;
              }
              k++;
            }
        }

        if(k!=32*size)
        {
          for(j=0;j<32*size;j++)
          {     
              if(j!=i-decal)
             {
                 int l;
                if(matrix[j][blocc]&( 1<<(31-restc)/*pow1(2,31-restc)*/))
                {
                  for(l=0;l<larg;l++)
                  {
                    matrix[j][l]= matrix[j][l] ^ matrix[i-decal][l];    
                  }
                  blocc2=(i-decal)/32;
                  restc2=(i-decal)%32;
  //                b[j/32]=b[j/32]^(!!(b[blocc2]&(1<<(31-restc2))/*pow1(2,31-restc2)*/)*(1<<(31-j%32)/*pow1(2,31-j%32)*/));
                }  
             }
          }
        }
 }
 /*ui kl=size*32-1;
 int solu=1;
 int count=0;
 ui *erreur;
 erreur=malloc(sizeof(ui)*size*32);
 int init;
 for(init=0;init<size*32;init++){erreur[init]=0;}
 while(!hamming_weight(matrix[kl],larg))
 {
     
      if(!!(b[kl/32]&(1<<(31-kl%32))) != 0)
      {
        int p;
        count++;
        for(p=0;p<size*32;p++)
        {
          if(init)
          {
            erreur[p]++;
          }
        }
        //print_vect(hist0[kl],size);  
      } 
      kl--;
 }

 int final=max(erreur,size);*/
return 0;
}

void* genere_with_10(ui **matrix,ui **sortie,ui *vect_i, ui *vect_o,int size)
{
  int i;
  for(i=0;i<32*size;i++)
  {
      int j;
      for(j=0;j<10;j++)
      {
        if(!!(i&pow1(2,j)))
        {
          vect_o[i/32 ]=vect_o[i/32]^( pow1(2,31-i%32) * !!( vect_i[j/32]&pow1(2,31-j%32)   ) );  
          xor_btb(sortie[i],matrix[j],sortie[i],SIZE_32);
        }
      }
  }
  return NULL;
}


void* gauss(ui **matrix, ui *b/*, ui *result*/,ui size)
{
//Gaussian elimination
 int i,j;
 ui blocc, restc;
 for(i=0;i<size*32;i++)
 {
        blocc=i/32;
        restc=i%32;
        if(!(matrix[i][blocc]&( pow1(2,31-restc))))
        {
            ui k=i+1;
            while(k<32*size)
            {
              if((matrix[k][blocc]&( pow1(2,31-restc))))
              { 
                  ui *tmp;
                  ui tmp1;
                  tmp=matrix[i];
                  matrix[i]=matrix[k];
                  matrix[k]=tmp;
                  tmp1=!!(b[blocc]&(pow1(2,31-restc))); 
                  //Can be improved
                  
                  b[blocc]=b[blocc]^(b[blocc]&pow1(2,31-restc))^((!!(b[k/32]&pow1(2,31-k%32)))*pow1(2,31-restc));
                  b[k/32]=(tmp1*pow1(2,31-k%32))^b[k/32]^(b[k/32]&(pow1(2,31- k%32)));
                  k=32*size;

              }
              k++;
            }
            if(k==32*size){printf("Fail"); return NULL;}
        }
        for(j=0;j<32*size;j++)
        {
             if(j!=i){
             int l;
            if(matrix[j][blocc]&( pow1(2,31-restc)))
             {
             for(l=0;l<size;l++)
             {
               matrix[j][l]= matrix[j][l] ^ matrix[i][l];    
             }
             b[j/32]=b[j/32]^ (!!(b[blocc]&pow1(2,31-restc))*(pow1(2,31-j%32)));
             }
             }
        }
 }
//Back substitution
/*result[SIZE_32-1]= b[SIZE_32-1]&1;
ui temp;
ui restic;
ui blocic;
for(i=SIZE_32*32-2;i>=0;i--)
{
  restic=i% 32;
  blocic=i/32;
  temp=!!(b[blocic]&pow1(2,31-restic));
  for(j=i+1;j<SIZE_32*32-1;j++)
  {
    restc = j % 32;
    blocc= j/ 32; 
    if(matrix[i][blocc]&( pow1(2,31-restc)))
    {temp +=  !!(result[blocc]&(pow1(2,31-restc)));}
  }
  result[blocic]= result[blocic]^((temp%2)*pow1(2,31-restic));
}*/
return NULL;
}

void* copy(ui **m1, ui**m2,int size1, int size2)
{
  int i,j;
    #pragma omp parallel for
  for(i=0;i<size1;i++)
  {
      for(j=0;j<size2;j++)
      {
        m2[i][j]=m1[i][j];  
      }
  }
}

void* copy_vect(ui *v1,ui *v2)
{
  int i;
  #pragma omp parallel for
  for(i=0;i<NB_IT;i++)
  {
        v2[i]=v1[i];  
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
       {v2[bloc]=v2[bloc]^(v2[bloc]&(1<<(31-rest)))^(!(v1[bloc]&pow1(2,31-rest))*(1<<(31-rest)));}
       else{v2[bloc]=v2[bloc]^(v2[bloc]&(1<<(31-rest)))^(v1[bloc]&pow1(2,31-rest));}

    }
} 

void* count_most_probable(ui ** mat, ui* b, float* sortie )
    //Warning sortie est de taille NB_IT 
{
  int i,j;
  ui w;
  for(i=0;i<NB_IT;i++)
  {
    ui bloc=i/32;
    ui rest=i%32;
    if(!!(b[bloc]&pow1(2,31-rest)))
    {
        w=hamming_weight(b,SIZE_32);
         
        for(j=0;j<NB_IT;j++)
        {
            ui blocj=j/32;
            ui restj=j%32;
            if(!!(mat[i][blocj]&pow1(2,31-restj)))
            {
                sortie[j]+= w;
            }
        }
    } 
  }
}

int geq(ui* a,ui* b)
{
  int res=1;
  int i;
  for(i=0;i<SIZE_32;i++)
  {
  res=res||!!(a[i]>=b[i]);
  }
  return res;
}



int leq(ui* a, ui*b)
{  int res=1;
    int i;
  for(i=0;i<SIZE_32;i++)
  {
  res=res||!!(a[i]<=b[i]);
  }
  return res;
}


void q_sort(ui **numbers, ui left, ui right)
{
  ui* pivot;
  ui  l_hold, r_hold;
 
  l_hold = left;
  r_hold = right;
  pivot = numbers[left];
  while (left < right)
      {
            while (geq(numbers[right], pivot) && (left < right))
                    right--;
              if (left != right)
                    {
                        //TODO : Possibility to use only a half:
                        xor_btb(numbers[left],numbers[left],numbers[left],SIZE_32);
                        xor_btb(numbers[left],numbers[right],numbers[left],SIZE_32);    
                                left++;
                                  }
                while (leq(numbers[left], pivot) && (left < right))
                        left++;
                  if (left != right)
                        {
                        //TODO : real affectation
                        xor_btb(numbers[right],numbers[right],numbers[right],SIZE_32);
                        xor_btb(numbers[right],numbers[left],numbers[right],SIZE_32);    
                                    right--;
                        }
                  }
  numbers[left] = pivot;
  ui pivot1;
  pivot1 = left;
  left = l_hold;
  right = r_hold;
  if (left < pivot1)
        q_sort(numbers, left, pivot1-1);
  if (right > pivot1)
        q_sort(numbers, pivot1+1, right);
}  


void quickSort(ui** numbers, ui array_size)
{
      q_sort(numbers, 0, array_size - 1);
}


ui search(ui** vect,ui pos)
{
  ui i=0;
  while(i<32*SIZE_32*3)
  {
    if(vect[pos][i/32]&(1<<(31-i%32)))
    {
      return i;
    }
   i++;
  }
  return (32*SIZE_32*3);
}

ui search2(ui** vect,ui pos)
{
  ui i=0;
  while(i<32*SIZE_32*3)
  {
    if((vect[pos+1][i/32]&(1<<(31-i%32)))&&(vect[pos][i/32]&(1<<(31-i%32))))
    {
      return i;
    }
   i++;
  }
  return (32*SIZE_32*3);
}


void* gauss_collide(ui** mat,ui* b, ui ** gen,ui pos)
{
 ui ** working;
 working=malloc(sizeof(ui*)*(32*SIZE_32+1));
 int in;
 for (in=0; in<32*SIZE_32+1; in++)
  {
     working[in] = malloc(sizeof(ui)*(SIZE_32+1));
  }
 zero_A(working,32*SIZE_32+1,SIZE_32+1);
 int init_id;
 for(init_id=0;init_id<SIZE_32*32+1;init_id++)
 {
   working[init_id][init_id/32]=working[init_id][init_id/32]^(pow1(2,31-init_id%32));   
 }
 ui i,j;
 ui blocc, restc;
 ui k=0;
 ui blocc2,restc2,blocc1,restc1;
 ui decal=0;
 for(i=0;i<SIZE_32*32;i++)
 {
       
        if(k==32*SIZE_32+1){decal++;}
        blocc=(i)/32;
        restc=(i)%32;
        k=i+1-decal;
        if(!(mat[i-decal][blocc]&(1<<(31-restc))))
        {
            while(k<32*SIZE_32+1)
            {
              if((mat[k][blocc]&(1<<(31-restc))))
              {
                  exchange(working,working,i-decal,k,SIZE_32+1); 
                  ui *tmp;
                  ui tmp1;
                  tmp=mat[i-decal];
                  mat[i-decal]=mat[k];
                  mat[k]=tmp;
                  blocc1=(i-decal)/32;
                  restc1=(i-decal)%32;
                  tmp1=!!(b[blocc1]&(1<<(31-restc1))); 
                  b[blocc1]=b[blocc1]^(b[blocc1]&(1<<(31-restc1)))^((!!(b[k/32]&(1<<(31-k%32))))*(1<<(31-restc1)));
                  b[k/32]=(tmp1*(1<<(31-k%32)))^b[k/32]^(b[k/32]&((1<<(31-k%32))));
                  k=32*SIZE_32+1;
              }
              k++;
            }
        }
        if(k!=32*SIZE_32+1)
        {
          for(j=0;j<32*SIZE_32+1;j++)
          {
             if(j!=i-decal)
             {
                int l;
                if(mat[j][blocc]&( 1<<(31-restc)))
                {
                  mix(working,working,j,i-decal,SIZE_32+1);
                  for(l=0;l<SIZE_32;l++)
                  {
                    mat[j][l]= mat[j][l] ^ mat[i-decal][l];    
                  }
                  blocc2=(i-decal)/32;
                  restc2=(i-decal)%32;
                  b[j/32]=b[j/32]^(!!(b[blocc2]&(1<<(31-restc2)))*(1<<(31-j%32)));
                }  
             }
          }
        }
 }
 zero_vect(gen[pos],3*SIZE_32); 
 
 xor_btb(gen[pos],working[SIZE_32*32-1],gen[pos],SIZE_32+1);
 //TODO : FREE WORKING!
 return NULL;
}


void* build3by2(ui** mat,ui* b,ui** gen,ui pos,ui* inj_1, ui* inj_2)
{
    ui pos_1,pos_2;     
    ui* temporary_delete;
    ui** work;
    ui* b1=malloc(sizeof(ui)*(32*SIZE_32 +1));
    work=malloc(sizeof(ui*)*(32*(SIZE_32)+1));
    int in;
     for (in=0; in<(1+(SIZE_32)*32); in++)
      {
         work[in] = malloc(sizeof(ui)*(SIZE_32));
      }
    for(in=0; in<(1+SIZE_32*32);in++)
    {
      b1[in]=b[in];
    }
    temporary_delete=malloc(sizeof(ui)*(SIZE_32));
    
    copy(mat,work, 32*SIZE_32, SIZE_32);
    zero_vect(work[32*SIZE_32],SIZE_32);
    xor_btb(work[32*SIZE_32],inj_1,work[32*SIZE_32], SIZE_32);    

    gauss_collide(work,b1,gen,pos);
   
    //TODO :mettre gen[pos]["SIZE_32+1"] dans gen[pos]["pos"]
   

    ui temp=!!(gen[pos][(32*SIZE_32 )/32]&(1<<(31-(32*SIZE_32)%32 )));
    gen[pos][(32*SIZE_32)/32]=gen[pos][(32*SIZE_32)/32]^(gen[pos][(32*SIZE_32)/32]&(1<<(31-(32*SIZE_32)%32 ))); 
    gen[pos][pos/32]=gen[pos][pos/32]^( temp*(1<<(31-pos%32)));

    pos_1=search(gen,pos); //find the first one in gen[pos] 
 
    copy(mat,work, 32*SIZE_32, SIZE_32);
    zero_vect(work[32*SIZE_32],SIZE_32);
    xor_btb(work[32*SIZE_32],inj_1,work[32*SIZE_32], SIZE_32);    
  
   
    zero_vect(temporary_delete,SIZE_32);
    xor_btb(work[pos_1],temporary_delete,temporary_delete, SIZE_32);    
   
    zero_vect(work[pos_1],SIZE_32);
    xor_btb(work[pos_1], inj_2 ,work[pos_1], SIZE_32);    
    for(in=0; in<(1+SIZE_32*32);in++)
    {
      b1[in]=b[in];
    }
    gauss_collide(work,b1,gen,pos+1);
    
    temp=!!(gen[pos+1][(32*SIZE_32)/32]&(1<<(31-(32*SIZE_32)%32 )));
    gen[pos+1][(32*SIZE_32)/32]=gen[pos+1][(32*SIZE_32)/32]^(gen[pos+1][(32*SIZE_32)/32]&(1<<(31-(32*SIZE_32)%32 ))); 
    gen[pos+1][pos/32]=gen[pos+1][pos/32]^( temp*(1<<(31-pos%32)));

    ui temp1=!!(gen[pos+1][(pos_1)/32]&(1<<(31-(pos_1)%32 )));
    gen[pos+1][(pos_1)/32]=gen[pos+1][(pos_1)/32]^(gen[pos+1][(pos_1)/32]&(1<<(31-(pos_1)%32 ))); 
    gen[pos+1][(pos+1)/32]=gen[pos+1][(pos+1)/32]^( temp*(1<<(31-(pos+1)%32)));


    pos_2=search2(gen,pos); //find the first one simultaneously in gen[pos] and gen[pos+1]
   //TODO : What if pos_2 doesn't exist? -> setjmp to retry? 
   //TODO : Use printf to see if there are lot of fails.     
   
    if (pos_2==32*3*SIZE_32){printf("\n\n\n\nARrrgl\n");}
    if(pos_2==pos_1){printf("\n\n\n bouh\n");}
    copy(mat,work, 32*SIZE_32, SIZE_32);
    zero_vect(work[32*SIZE_32],SIZE_32);
    xor_btb(work[32*SIZE_32],inj_1,work[32*SIZE_32], SIZE_32);    
  

    zero_vect(work[pos_2],SIZE_32); 
    xor_btb(work[pos_2],temporary_delete ,work[pos_2], SIZE_32);   
    for(in=0; in<(1+SIZE_32*32);in++)
    {
      b1[in]=b[in];
    }
    gauss_collide(work,b1,gen,pos+2);
  
    temp=!!(gen[pos+2][(SIZE_32*32 )/32]&(1<<(31-(32*SIZE_32)%32 )));
    gen[pos+2][(SIZE_32*32)/32]=gen[pos+2][(32*SIZE_32)/32]^(gen[pos+2][(32*SIZE_32)/32]&(1<<(31-(32*SIZE_32)%32 ))); 
    gen[pos+2][pos/32]=gen[pos+2][pos/32]^( temp*(1<<(31-pos%32)));

    temp1=!!(gen[pos+2][(pos_1)/32]&(1<<(31-(pos_1)%32 )));
    gen[pos+2][(pos_1)/32]=gen[pos+2][(pos_1)/32]^(gen[pos+2][(pos_1)/32]&(1<<(31-(pos_1)%32 ))); 
    gen[pos+2][(pos+1)/32]=gen[pos+2][(pos+1)/32]^( temp1*(1<<(31-(pos+1)%32)));


    ui temp2=!!(gen[pos+2][(pos_2)/32]&(1<<(31-(pos_2)%32 )));
    gen[pos+2][(pos_2)/32]=gen[pos+2][(pos_2)/32]^(gen[pos+2][(pos_2)/32]&(1<<(31-(pos_2)%32 ))); 
    gen[pos+2][(pos+2)/32]=gen[pos+2][(pos+2)/32]^( temp2*(1<<(31-(pos+2)%32)));

//FREE WORK
    return NULL;
}




void* build_eqs(ui** mat,ui* b)
{
  ui** gen ;
  ui i;
  ui* inj_1;
  ui* inj_2;
  gen= malloc(sizeof(ui*)*32*3*(SIZE_32));
  inj_1=malloc(sizeof(ui)*SIZE_32);
  inj_2=malloc(sizeof(ui)*SIZE_32);
  int in;
     for (in=0; in<(32*3*SIZE_32); in++)
      {
         gen[in] = malloc(sizeof(ui)*3*(SIZE_32));
      }
  zero_A(gen,32*3*SIZE_32,3*SIZE_32);

  for(i=0;i<SIZE_32*32;i++)
  {
    zero_vect(inj_1,SIZE_32);
    xor_btb(mat[2*i+32*SIZE_32],inj_1,inj_1,SIZE_32);
    zero_vect(inj_2,SIZE_32);
    xor_btb(mat[2*i+1+32*SIZE_32],inj_2,inj_2,SIZE_32);

    build3by2(mat,b,gen,3*i,inj_1,inj_2);
  }
  gauss2(gen,b,3*SIZE_32,3*SIZE_32);
  print_mat(gen,32*3*SIZE_32,3*SIZE_32);   
}



unsigned long binom(ui n, ui p)
{
  unsigned long binom;
  unsigned long k;			
  if (p > n-p) 
    p = n-p;
  for (binom = 1, k = 1; k <= p; k++)
    binom = (binom * (n-p+k)) / k;
  return binom;
}

void* generate_hamming(ui** work,ui ndepart,ui n,ui w,unsigned long begin, unsigned long end)
/* Generate all n-words of hamming weight w in work.
 * begin and end are initialized with 0 and binomial(w,n) respectively
 * ndepart=n initially. Actually ndepart is a constant
 */  
{
    unsigned long pos;
    if(w==0||n==0){return NULL;}
    if(n==w)
    {
       unsigned long j,i;
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
        unsigned long i;
        for(i=begin;i<begin+pos;i++)    
        {
            work[i][(ndepart-n)/32]=work[i][(ndepart-n)/32]^(1<<(31-(ndepart-n)%32));
        }
        generate_hamming(work,ndepart,n-1,w-1,begin,begin+pos);
        generate_hamming(work,ndepart,n-1,w,begin+pos,end);
    }
}

int geq1(ui* a,ui* b,ui begin,ui end)
/*(a) Greater or equal (b), considering only
 * the bits between begin and end. 
 * It doesn't work for all bits, maximum allbits-1
 * but pratically it's okay*/
{
  int res=1;
  int i;
  for(i=0;i<SIZE_32;i++)
  {
      if(begin/32<i)
      {
          if(end/32>i)
              res=res&&!!(a[i]>=b[i]);
          else if(end/32==i)
              res=res&&!!((a[i]&(((1<<(end%32))-1)<<(31-end)))>=(b[i]&(((1<<(end%32))-1)<<(31-end))));
      }
      else if(begin/32==i)
      {
          if(end/32>i)
          {
              res=res&&!!((a[i]&(((1<<(-(begin%32)%32))-1)))>=(b[i]&(((1<<(-(begin%32)%32))-1))));
          }
          else if(end/32==i)
          {
              res=res&&!!((a[i]&(((1<<((end-begin)%32))-1)<<(31-end+1)))>=(b[i]&(((1<<((end-begin)%32))-1)<<(31-end+1))));
          }
      }
 }
  return res;
}


int leq1(ui* a,ui* b,ui begin,ui end)
{
  int res=1;
  int i;
  for(i=0;i<SIZE_32;i++)
  {
      if(begin/32<i)
      {
          if(end/32>i)
              res=res&&!!(a[i]<=b[i]);
          else if(end/32==i)
              res=res&&!!((a[i]&(((1<<(end%32))-1)<<(31-end)))<=(b[i]&(((1<<(end%32))-1)<<(31-end))));
      }
      else if(begin/32==i)
      {
          if(end/32>i)
              res=res&&!!((a[i]&(((1<<(-(begin%32)%32))-1)))<=(b[i]&(((1<<(-(begin%32)%32))-1))));
          else if(end/32==i){
              res=res&&!!((a[i]&(((1<<((end-begin)%32))-1)<<(31-end+1)))<=(b[i]&((((1<<((end-begin)%32))-1)<<(31-end+1)))));}
      }
 }
  return res;
}


void echanger(ui** tableau,ui a, ui b)
{
        ui* temp = malloc(sizeof(ui)*SIZE_32);
        ui* zero = malloc(sizeof(ui)*SIZE_32);
        zero_vect(zero,SIZE_32);
        zero_vect(temp,SIZE_32);
        xor_btb(zero,tableau[a],temp,SIZE_32);
        xor_btb(zero,tableau[b],tableau[a],SIZE_32);
        xor_btb(zero,temp,tableau[b],SIZE_32);
        free(temp);
        free(zero);
}
 
void quickS(ui** tableau,ui* b, long debut, long fin,ui begin,ui end)
{
      long gauche = debut-1;
      long droite = fin+1;
      ui* pivot = tableau[debut];
      if(debut >= fin)
                  return;
       
          while(1)
          {
           do {droite--;}while(!geq1(tableau[droite], pivot,begin,end));
           
           do gauche++;while(!leq1(tableau[gauche], pivot,begin,end));
           if(gauche < droite)
           {
                  ui swap;
                  swap=!!(b[gauche/32]&(1<<(31-gauche%32)));
                  b[gauche/32]=b[gauche/32]^(b[gauche/32]&(1<<(31-gauche%32)))^
                              ((!!(b[droite/32]&(1<<(31-droite%32))))*(1<<(31-gauche%32)));
                  b[droite/32]=(swap*(1<<(31-droite%32)))^b[droite/32]^(b[droite/32]&((1<<(31-droite%32))/*pow1(2,31- k%32)*/));
                  echanger(tableau, gauche, droite);
           }
           else break;
          }
              quickS(tableau,b, debut, droite,begin,end);
              quickS(tableau,b, droite+1, fin,begin,end);
}

void quicky(ui** numbers,ui* b, long array_size,ui begin,ui end)
{
      quickS(numbers,b, 0, array_size - 1,begin,end);
}

int equal(ui* a,ui* b)
{
  int res=1,i;
  for(i=0;i<SIZE_32;i++)
      res=res&&(a[i]==b[i]);
  return res;
          

}

void* BKW_reduction(ui** work,ui* snd,ui a, ui b,ui sizel)
/*This function make (a) steps of size (b)
 * to reduce the size of work. The noise increase
 * exponentially
 */
{
  ui ba;
  int i;
  int asuppr=-1;
  ui* zero = malloc(sizeof(ui)*SIZE_32);
  ui* temp=malloc(sizeof(ui)*SIZE_32);
  zero_vect(temp,SIZE_32);
  zero_vect(zero,SIZE_32);
  for(ba=0;ba<a;ba++)
  {
      quicky(work,snd,sizel,ba*b,(ba+1)*b);
      for(i=0;i<sizel;i++)
      {
          if(!(equal(work[i],zero)))
          {
          if(leq1(temp,work[i],ba*b,(ba+1)*b)&&geq1(temp,work[i],ba*b,(ba+1)*b))
          {
              snd[i/32]=snd[i/32]^((!!(snd[asuppr/32]&(1<<(31-asuppr%32))))*(1<<(31-i%32 )));
              xor_btb(work[i],temp,work[i],SIZE_32);
          }
          else
          {
              if(asuppr>=0)zero_vect(work[asuppr],SIZE_32);
              asuppr=i;
              xor_btb(work[i],zero,temp,SIZE_32);
          }
          }
      }    
  }
  ui from_end=1;
  for(i=0;i<CSTVAUDENAY;i++)
  {
      if(equal(work[i],zero))
        {
            //TODO: BEURK : DIVIDE INSTRUCTION BY 2 : IT'S JUST
            //A PERMUTATION OF (FROM_END) AND (I) !!!
            zero_vect(temp, SIZE_32);
            xor_btb(zero,work[i],temp,SIZE_32);
            zero_vect(work[i],SIZE_32);
            xor_btb(zero,work[sizel-from_end],work[i],SIZE_32);
            zero_vect(work[sizel-from_end],SIZE_32);
            xor_btb(zero,temp,work[sizel-from_end],SIZE_32);
            ui blocc1,restc1,tmp1;
            blocc1=(i)/32;
            restc1=(i)%32;
            tmp1=!!(snd[blocc1]&(1<<(31-restc1))); 
            snd[blocc1]=snd[blocc1]^(snd[blocc1]&(1<<(31-restc1)))^((!!(snd[(sizel-from_end)/32]&(1<<(31-(sizel-from_end)%32))))*(1<<(31-restc1)));
            snd[(sizel-from_end)/32]=(tmp1*(1<<(31-(sizel-from_end)%32)))^snd[(sizel-from_end)/32]^(snd[(sizel-from_end)/32]&((1<<(31-(sizel-from_end)%32))));
            from_end++;
        }
  }
  free(zero);
  free(temp);
  return NULL;
}


float AKMV_estimator(ui** bici,ui* faibi,ui* faici,ui* c,float twotothe2k)
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
        if ((times%2))
        {
        result+= !!(faibi[n/32]&(1<<(31-n%32)))
                *!!(faici[n/32]&(1<<(31-n%32)))
                *(-1);
        }
        else
        {
        result+= !!(faibi[n/32]&(1<<(31-n%32)))
                *!!(faici[n/32]&(1<<(31-n%32)));
        } //CECI EST DEPOURVU DE SENS
    }
    result=result*twotothe2k/CSTVAUDENAY;
    return result;
}

void* generate_4_AKMV(ui** work,ui* b,ui** bici,ui* faibi,ui* faici)
/*This function makes the datas for akmv estimator*/
{
    //TODO: put 0-lines at the end of (work)
    int n;
    for(n=0;n<CSTVAUDENAY;n++)
    {
        faici[n/32]=faici[n/32]^((1<<(31- n%32 ))*!!(b[(2*n)/32]&(1<<(31-(2*n)%32))));
        faibi[n/32]=faibi[n/32]^((1<<(31- n%32 ))*!!(b[(2*n+1)/32]&(1<<(31-(2*n+1)%32))));
        xor_btb(work[2*n+1],work[2*n],bici[n],SIZE_32);
    }
}

void* bruteforce_attack(ui** hamming,ui* b,ui** samples,ui sizebrute)
{
    int i;
    float value;
    ui* faibi = malloc(sizeof(ui)*CSTVAUDENAY);
    ui* faici = malloc(sizeof(ui)*CSTVAUDENAY);
    ui* bici1 = malloc(sizeof(ui)*CSTVAUDENAY*SIZE_32);
    long* bici=malloc(sizeof(long)*CSTVAUDENAY);
    for(i=0;i<CSTVAUDENAY;i++)
      bici[i]=(long)&(bici1[SIZE_32*i]);  
    generate_4_AKMV(samples, b,(ui**)bici,faibi,faici);
    float local;
    local=0;
   //TODO : Make all estimators in one shot? 
   #pragma omp parallel for 
    for(i=0;i<sizebrute;i++)
    {
      value = AKMV_estimator((ui**)bici,faibi,faici,hamming[i],1<<30);
      if(value>local){printf("%f\n",value);print_vect(hamming[i],SIZE_32);local=value;}
    }
    printf("%f",value);
    free(faici);
    free(faibi);
    free(bici); 
    free(bici1);
}

int main (void)
{
      //  //omp_set_num_threads(4);
        srand(time(NULL));
      //  int size=4*SIZE_32;
      //  int nbrpermut=2;
      //  ui **A;
      //  ui **A2;
      //  ui **A3;
      //  ui **working;
      //  ui **working2;
      //  ui **working4;
      //  ui **working3;
      //  A2 = malloc(sizeof(ui*)*32* size); //WARNING MODIFICATION
      //  int i;
      //  for (i=0; i<32*size; i++)
      //  {
      //        A2[i] = malloc(sizeof(ui)*SIZE_32);
      //  }
      //  A3 = malloc(sizeof(ui*)*32* size); //WARNING MODIFICATION
      //  for (i=0; i<32*size; i++)
      //  {
      //        A3[i] = malloc(sizeof(ui)*SIZE_32);
      //  }
      //  A = malloc(sizeof(ui*)*32* size); //WARNING MODIFICATION
      //  for (i=0; i<32*size; i++)
      //  {
      //        A[i] = malloc(sizeof(ui)*SIZE_32);
      //  }
      //  working = malloc(sizeof(ui*)*32* size); //AND HERE
      //  for (i=0; i<32*size; i++)
      //  {
      //     working[i] = malloc(sizeof(ui)*size);
      //  }
      //  working2 = malloc(sizeof(ui*)*32* size); //AND HERE
      //  for (i=0; i<32*size; i++)
      //  {
      //     working2[i] = malloc(sizeof(ui)*size);
      //  }
      //  working4 = malloc(sizeof(ui*)*32* size); //AND HERE
      //  for (i=0; i<32*size; i++)
      //  {
      //     working4[i] = malloc(sizeof(ui)*size);
      //  }
      //  working3 = malloc(sizeof(ui*)*32*nbrpermut*(size-SIZE_32)); //AND HERE
      //  for (i=0; i<32*nbrpermut*(size-SIZE_32); i++)
      //  {
      //     working3[i] = malloc(sizeof(ui)*size);
      //  }
      //  unsigned int *b;
      //  unsigned int *result;
      //  unsigned int *result1;
      //  unsigned int *result2;
      //  ui *result3;
      //  ui *result4;
      //  /*Initializationn*/
      //  b=malloc(sizeof(ui)*SIZE_32);
      //  result1=malloc(sizeof(ui)*size);
      //  result2=malloc(sizeof(ui)*size);
      //  result4=malloc(sizeof(ui)*size);
      //  result3=malloc(sizeof(ui)*nbrpermut*(size-SIZE_32));
      //  result=malloc(sizeof(ui)*size);
      //      //  /*End initialization */
      //      //  
      //      // //print_vect(result,size);
      //      //  //print_vect(result1,size); 
      //      //     //print_vect(result,size);
      //      //     //zero_A(working);
      //      //     //genere_with_10(A,working,result,result2);
      //      //     //times(working,b,result1,SIZE_32); 
      //      //     //print_vect(result,SIZE_32);
      //      //     //print_vect(result1,SIZE_32);
      //      //   int gagne=0;
      //      // //while(!gagne||hamming_weight(result3,size)>60)
      //      // //{
      //      //             zero_A(working,size*32,size);
      //      //  int init_id;
      //      //  for(init_id=0;init_id<size*32;init_id++)
      //      //  {
      //      //    working[init_id][init_id/32]=working[init_id][init_id/32]^(pow1(2,31-init_id%32));   
      //      //  }
      //      // 
      //      //  generate_A(A,size*32,SIZE_32);
      //      //  generate_vect(b,SIZE_32);
      //      //  zero_vect(result,size);
      //      //  times(A,b,result,size); 
      //      //  add_noise(result,result1,0.05,size);
      //      //  print_vect(result,size);
      //      //  print_vect(result1,size);
      //      // // add_noise(result1,result2,0.0,size);
      //      //  
      //      //  //copy(working,working2,size*32,size);
      //      //  //randomize(A2,result2,working2,size);
      //      //  //randomize(A,result1,working,size);
      //      //  //gagne=gauss_test(A2,result2,working2,size);
      //      //  //print_mat(working2,size*32,size);
      //      //  //gagne=gauss_test(A,result1,working,size);
      //      //  int bouc;
      //      //  for(bouc=0;bouc<nbrpermut;bouc++)
      //      //  {
      //      //  add_noise(result1,result4,0.0,size);
      //      //  copy(A,A3,size*32,SIZE_32);
      //      //  copy(working,working4,size*32,size);
      //      //  randomize(A3,result4,working4,size);
      //      //  //zero_vect(result3,nbrpermut*(size-SIZE_32));
        //      //  gagne=gauss_test(A3,result4,working4,size);
      //      //  for(gagne=0;gagne<32*(size-SIZE_32);gagne++)
      //      //    {
      //      //      xor_btb(working3[bouc*32*(size-SIZE_32)+gagne],
      //      //              working3[bouc*32*(size-SIZE_32)+gagne],working3[bouc*32*(size-SIZE_32)+gagne],size);
      //      //      xor_btb(working3[bouc*32*(size-SIZE_32)+gagne],working4[SIZE_32*32+gagne],working3[bouc*32*(size-SIZE_32)+gagne],size);
      //      //      result3[(bouc*32*(size-SIZE_32)+gagne)/32]=result4[SIZE_32+(gagne)/32];
      //      //    }
      //      //  }
      //      //  gauss2(working3,result3,nbrpermut*(size-SIZE_32),size);
      //      //  print_mat(working3,32*(size-SIZE_32)+3,size);
      //      //  print_vect(result3,size);
      //      //  fflush(stdout);
      //      ////}
      //      //  gagne=1;
      //      //  //print_vect(result3,3*(size-SIZE_32)); 
      //      //  printf("\n %d \n",gagne); 
      //      //  return 0;
      //      //  // print_vect(result2,size);
      //      //     //print_vect(result1,SIZE_32);
      //      //     // exhaustive(A,b,result);
      //  generate_A(A,size*32,SIZE_32);
      //  generate_vect(b,SIZE_32);
      //  zero_vect(result,size);
      //  times(A,b,result,size); 
      //  add_noise(result,result1,0.05,size);
      //  build_eqs(A,result1);

  int i;
  int nberreur=3;
  unsigned long t=binom(32*SIZE_32,nberreur);
  long* work =malloc(sizeof(long)*NB_IT);
  ui* work1= malloc(sizeof(ui)*NB_IT*(SIZE_32));
  long* hamming =malloc(sizeof(long)*t);
  ui* hamming1= malloc(sizeof(ui)*SIZE_32*t);
  ui* s = malloc(sizeof(ui)*SIZE_32);
  ui* b = malloc(sizeof(ui)*NB_IT/32); 
  for(i=0;i<t;i++)
  {
     hamming[i]=(long)&(hamming1[SIZE_32*i]); 
  }
  for(i=0;i<NB_IT;i++)
  {
     work[i]=(long)&(work1[SIZE_32*i]); 
  }
  generate_A((ui**)work,NB_IT,SIZE_32);
  zero_vect(b,NB_IT/32);
  s[0]=(1<<31)+(1<<14)+(1<<1);
  times((ui**)work,s,b,NB_IT/32);
  add_noise(b,b,0.3,NB_IT/32);
  print_vect(s,SIZE_32);
  generate_hamming((ui**)hamming,SIZE_32*32,SIZE_32*32,nberreur,0,t);
 // BKW_reduction((ui**)work,(ui*)b,5,21,NB_IT);
  bruteforce_attack((ui**)hamming,(ui*)b,(ui**)work, t);
  free(work1);
  free(work);
  free(hamming1);
  free(hamming);
  free(s);
  free(b);
  return 0;	
}
