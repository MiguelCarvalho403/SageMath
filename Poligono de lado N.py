n = 3
alfa = 2*pi/n
l = 5
r = l/(sqrt(2 - 2* cos(alfa) ))
p_c = []

for i in range(1, n + 1):
    m = pi/2 + alfa * (i - 1)
    print(m)
    p_t = [r*cos(m), r*sin(m)]    
    p_c.append(p_t)

print(p_c)
print(3*cos(90)/sqrt(-2*cos(90) + 2))
p = polygon(p_c, color='red')
p.show(xmin = -r, xmax = r, ymin = -r, ymax = r)
