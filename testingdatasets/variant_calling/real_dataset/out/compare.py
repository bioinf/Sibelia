gdvar = [int(line.strip().split()[0]) for line in open('truevariant.txt')]
sibvar = [int(line.strip().split()[0]) + 1 for line in open('variant.txt')]

gdhit = []
sibhit = []
nearhit = []
for v in sibvar:
	if v in gdvar:
		gdhit.append(v)
		sibhit.append(v)
	else:
		for m in gdvar:
			if abs(m - v) < 10:
				gdhit.append(m)
				nearhit.append(v)
				break
for v in sibhit:
	print 'Hit!', v

for v in nearhit:
	print 'Near hit!', v

for v in gdvar:
	if v not in gdhit:
		print 'Not hit!', v

for v in sibvar:
	if v not in sibhit and v not in nearhit:
		print 'False hit!', v