import sys
from matplotlib import pyplot as plt
for e in range(len(sys.argv)):
    if sys.argv[e] == "-i":
        pos1 = e + 1

file = open(sys.argv[pos1] + '.fa','r')
protein = file.read()
a = protein.find('\n')
b= protein.find('>',a)
c= protein.find('\n',b)
protein_1 = protein[a+1:b-2]
protein_2 = protein[c+1:len(protein)-1]
protein_1.replace("\n","")
protein_2.replace("\n","")
np =0
#protein_1 = "GATC"
#protein_2 = "GAGC"
one = len(protein_1)    
two = len(protein_2)
one_copy = one 
two_copy = two	
def fill_0 (one ,two , m):
	for i in range(0,two+1):
		matrix.append([])
		for j in range(0,one+1):
			matrix[i].append(0)

def match_score (protein_1,protein_2,i,j,ans,match,mismatch):
		
		if protein_1[j-1] == protein_2[i-1]:
			ans = match 

		else :
			ans = mismatch		

matrix = []

fill_0(one,two,matrix)

#print(matrix)
one_allign = ""
two_allign = ""
#scores for match gap and mismatch
gap_penalty = 0
match =1
mismatch = -1

#adding gap penalty in first most row and column
for j in range(1,one+1):
	matrix[j][0] = gap_penalty *j 
for i in range(1,two+1):
	matrix[0][i] = gap_penalty *i


#print(matrix)#######################################

# filling matrix 
diagonal_score =0 

for i in range(1,two+1):
	for j in range(1,one+1):
		gap_one = matrix[i][j-1]
		gap_one += gap_penalty
		gap_two = matrix[i-1][j] 
		gap_one += gap_penalty
		#match_score(protein_1,protein_2,i,j,diagonal_score,match,mismatch)
		if protein_1[j-1] == protein_2[i-1]:
			diagonal_score = match 

		else :
			diagonal_score = mismatch
		
		#print(diagonal_score)
		
		diagonal_match = matrix[i-1][j-1] +diagonal_score 
		ans = max(gap_one,gap_two,diagonal_match)

		matrix[i][j] = ans

# backtracing


while one_copy >0 or two_copy >0 :
	#print(one_copy,two_copy)
	if protein_1[one_copy-1 ] == protein_2[two_copy -1]:
		np = match
	else :
		np = mismatch
		pass


	if  matrix[two_copy-1][one_copy-1]== matrix[two_copy][one_copy] -np :
		#if ( one_copy-1 >0):
		one_allign += protein_1[one_copy-1]
		one_copy -=1

		if ( one_copy-1 <= 0):
			two_allign+= "_"
			two_copy-=1
		else:	
			two_allign += protein_2[two_copy-1]
			two_copy -=1

	elif matrix[two_copy][one_copy-1]== matrix[two_copy][one_copy] - gap_penalty:
		
		one_allign += protein_1[one_copy-1]
		one_copy -=1

		two_allign += "_"
		
	elif matrix[two_copy-1][one_copy] == matrix[two_copy][one_copy] - gap_penalty:
		one_allign += "_"
		if ( one_copy-1 <= 0):
			two_allign+= "_"
			two_copy-=1
		else:	
			two_allign += protein_2[two_copy-1]
			two_copy -=1
	else : 
		pass	

sum_matrix = matrix
one_allign = one_allign [::-1]
two_allign = two_allign[::-1]


print("one_align")
print(one_allign)
print("two_align")
print(two_allign)
print("sum_matrix")
print(sum_matrix)		

print("dotplot")
print(protein_1)   #columns
print(protein_2)   #rows
# protein_1 = "abcd"
# protein_2 = "dbct"
columns = len(protein_1)
rows = len(protein_2)
dotplot = [[0]*columns for i in range(rows)]

for i in range (rows):
	for j in range (columns):
		if(protein_2[i] == protein_1[j]):
			dotplot[i][j] = 1
		else:
			dotplot[i][j] = 0


for i in range (rows):
	print(dotplot[i])
plt.imshow(dotplot)
plt.show()

file.close()