n=64;
nline=3*n;
epsilon=0.048;
s=random_matrix(GF(2),n,1);
a=random_matrix(GF(2),nline,n);
nu=random_matrix(GF(2),nline,1,density=epsilon)
b=a*s+nu
i=0
for i in range(nline):
  if b[i,:]==1:
    j=0
    while (j < n) and (s[j]==0 or a[i,j]==0):
      j += 1
    a[i,j]=0
b=a*s+nu
breal=a*s
#def take_2_build_3(mat,a,b,i,j,l):
#  result = matrix(GF(2),0,(l+1)*n)
#  pump1 = block_matrix([ mat , a ],ncols=1)
#  vect1 = pump1.kernel().basis()[0] #Kernel of dimension 1
#  result = block_matrix([ result , block_matrix( [ matrix(vect1),matrix(GF(2),1,l*n-1) ],ncols=2) ],ncols=1)
#  if i!=n:
#    result[0, i]=result[0, n]
#    result[0, n]=0
#  compt=0
#  while(vect1[compt]==0):
#    compt+=1
#  pump2 = block_matrix([mat[0:compt,:],b, mat[compt+1:n],a],ncols=1)
#  vect2 = pump2.kernel().basis()[0] #Kernel of dimension 1
#  if vect2[n]==0:
#    raise Exception('Try again') #Proba 1/2 -> origin of fail
#  result = block_matrix([ result , block_matrix( [ matrix(vect2),matrix(GF(2),1,l*n-1) ],ncols=2) ],ncols=1)
#  result[1,j]=result[1,compt]  
#  result[1,compt]=0
#  if i!=n:
#    result[1, i]=result[1, n]
#    result[1, n]=0
#  pump3 = block_matrix([mat[0:compt,:],b,mat[compt+1:n],mat[compt,:]],ncols=1)
#  vect3 = pump3.kernel().basis()[0]  #Kernel of dimension 1
#  result = block_matrix([ result , block_matrix( [ matrix(vect3),matrix(GF(2),1,l*n-1) ],ncols=2) ],ncols=1)
#  result[2, j]=result[2, compt]
#  result[2, compt]=result[2, n] 
#  result[2, n]=0
#  return result
    
#def make_mat(mat,var,p):
#  result=matrix(GF(2),n,(p+1)*n)
#  i=0
#  j=0
#  while i<var:
#    try:
#      result=block_matrix([result,take_2_build_3(mat[0:n,:],mat[n+j+2*i,:],mat[n+j+2*i+1,:],n+2*i,n+2*i+1,p)],ncols=1)
#      i += 1
#    except Exception as ins:
#      j += 1
#  return result
     

def make_and_reduce(a,b):
  noyau=a.kernel()
  res=(noyau.matrix())
  reduit=(res.transpose().kernel().matrix().change_ring(ZZ).LLL())
  return reduit.change_ring(GF(2)).change_ring(ZZ)[0,:]*matrix.ones(3*n,1)
  
