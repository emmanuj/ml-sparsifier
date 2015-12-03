import math

def spearman_rank(X, Y):
	if(len(X) != len(Y)):
		raise ValueError("The two list have to be equal")

	x = rank(X)
	y = rank(Y)

	d = 0.0
	for i,j in zip(x,y):
		d += math.pow((i-j),2)

	n = len(X)

	return 1 - ((6 * d)/(n * (math.pow(n,2) - 1)))

def rank(data): # 4.5, 2.1, 6.7
	results = []
	sorted_d = sorted(data) # 2.1,4.5,6.7
	ranking = [i for i in range(1, len(data)+1)]
	#print sorted_d, ranking
	#take care of ties
	i=0
	while i < len(sorted_d):
		#print ranking
		j = i + 1
		sum_d= ranking[i]
		count = 1
		while j < len(sorted_d) and (sorted_d[i] == sorted_d[j]):
			sum_d+=ranking[j]
			count+=1
			j+=1

		k = i
		while k < (i + count):
			ranking[k] = sum_d/float(count)
			k+=1
		#print ranking

		i+=count
	#print ranking
	dkey = dict(zip(sorted_d, ranking))
	#print ranking
	#print dkey
	for d in data:
		results.append(dkey[d])
	return results
'''
def test():
	# test rank
	l = [56,75,45,71,61,64,58,80,76,61]
	l2 = [66,70,40,60,65,56,59,77,67,63]

	print spearman_rank(l,l2)

if __name__ == '__main__':
	test()
'''
