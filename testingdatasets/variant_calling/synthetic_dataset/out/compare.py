gdline = [line for line in open('truevariant.txt')]
gdvar = [int(line.strip().split()[0]) for line in gdline]
sibvar = [int(line.strip().split()[0]) + 1 for line in open('variant.txt')]
varlen = dict()
lenhit = dict()
for line in gdline:
	line = line.strip().split()
	length = max(len(line[1]), len(line[2]))
	varlen[int(line[0])] = length
	lenhit[length] = 0



gdhit = []
sibhit = []
nearhit = []
for v in sibvar:
	if v in gdvar:
		gdhit.append(v)
		sibhit.append(v)
	else:
		for m in gdvar:
			if abs(m - v) < 20:
				gdhit.append(m)
				sibhit.append(v)
				#nearhit.append(v)
				break
for v in gdhit:
	print 'Hit!', v
	vlen = varlen[v]
	lenhit[vlen] = lenhit[vlen] + 1

for v in nearhit:
	print 'Near hit!', v

for v in gdvar:
	if v not in gdhit:
		print 'Not hit!', v

for v in sibvar:
	if v not in sibhit and v not in nearhit:
		print 'False hit!', v

print 'Exact match length hit:'
print lenhit