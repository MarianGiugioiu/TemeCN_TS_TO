f = open("Input", "r")
x = f.read()
y = ""
for c in x:
    s = c
    if c == 'ă':
        print(c)
        s = "a"
    y += s
#print(y)