

txt = '''Li: -10.2676 3
B : -10.0506 5
C : -3.5229  6
N : -4.1549  7
O : -3.3979  8
Ne: -4.2218  10
Na: -6.5229  11
Mg: -5.5229  12
Al: -6.6990  13
Si: -5.3979  14
P : -6.7959  19
S : -5.0000  16
Cl: -7.000   17
Ar: -5.5229  18
K : -7.9586  19
Ca: -7.6990  20
Ti: -9.2366  22
V : -10.0000 23
Cr: -8.0000  24
Mn: -7.6383  25
Fe: -5.5229  26
Ni: -7.0000  28
Cu: -8.8239  29
Zn: -7.6990  30'''
#txt = txt.replace("\n"," ")
#abuns = filter(None,txt.split(" "))
#print abuns
X = 1.0
Y = 10.0**-1.0223
Z = 0.0
lines = txt.split("\n")
for line in lines:
    line = line.replace(":","")
    items = line.split(" ")
    items = filter(None,items)
    ratio = float(items[1])
    mass = float(items[2])
    print ratio, mass, (10.0**ratio)*mass
    Z += (10.0**ratio)*mass
print Z