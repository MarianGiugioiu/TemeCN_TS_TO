f = open("InputLPDA", "r")
Q = f.readline().split() #starile
M = f.readline().split() #alfabetul de intrare
P = f.readline().split() #alfabetul de iesire
Z = f.readline().split() #alfabetul stivei
q0 = f.readline()[:-1] #starea initiala
z0 = f.readline()[:-1] #simbolul initial al stivei
F = f.readline().split() #starile finale
N = int(f.readline()) #numar tranzitii
D = {} #tranzitiile
for i in range(N):
    x = f.readline().split()
    #se va folosi simbolul . pentru lambda
    a = tuple(x[:3]) #stare initiala, simbol imput, simbol de scos din stiva
    b = tuple(x[3:]) #stare finala, simbol output, simboluri de adaugat in stiva
    if not (a in D):
        D[a] = []
    D[a].append(b)
    #print(D[a])

def show_path(str, path):
    print('(',q0,str[:-1],z0, end =" ) |- ")
    for t in path:
        (a,b,k,S) = t
        print('(', b[0],str[k:-1],''.join(S)[::-1], end =" ) |- ")
    print("accept")

def solve(str, k, q, S, sol, path):
    if k == len(str) - 1 and q in F: #daca am terminat cuvantul si sunt in stare finala
        print(sol)
        #show_path(str,path) #afiseaza drumul
    else:
        if S: #mai sunt simboluri in stiva
            a = (q, str[k], S[-1]) #aleg tranzitiile care au nevoie de litera curenta din cuvant si scot simbolul din varful stivei
            #print(a)
            #print(S)
            #print("##")
            if a in D:
                for b in D[a]:
                    #print(b)
                    S_current = S.copy()
                    sol_current = sol
                    path_current = path.copy()

                    del S_current[-1] #scot simbolul din varful stivei
                    if(b[2] != '.'): #adaug daca este nevoie simboluri in stiva
                        S_current += list(b[2][::-1])
                    if(b[1] != '.'): #adaug daca este nevoie simbol la output
                        sol_current += b[1]
                    #daca parcurg litera curenta in cuvant trebuie sa cresc k
                    if k < len(str) - 1:
                        path_current.append((a, b, k+1, S_current)) #adaug tranzitia in path
                        solve(str, k+1, b[0], S_current, sol_current,path_current)

            a = (q, '.', S[-1]) #aleg tranzitiile care au lambda in loc de litera curenta din cuvant si scot simbolul din varful stivei
            # print(D[a])
            if a in D:
                for b in D[a]:
                    #print(b)
                    S_current = S.copy()
                    sol_current = sol
                    path_current = path.copy()
                    del S_current[-1]
                    if (b[2] != '.'):
                        S_current += list(b[2][::-1])
                    if (b[1] != '.'):
                        sol_current += b[1]
                    #aici nu trebuie sa cresc k
                    path_current.append((a, b, k, S_current))
                    solve(str, k, b[0], S_current, sol_current, path_current)

        a = (q, str[k], '.') #aleg tranzitiile care au nevoie de litera curenta din cuvant si nu scot simbolul din varful stivei
        # print(D[a])
        if a in D:
            for b in D[a]:
                #print(b)
                S_current = S.copy()
                sol_current = sol
                path_current = path.copy()
                #del S_current[-1] #nu mai este nevoie sa scot varful listei
                if (b[2] != '.'):
                    S_current += list(b[2][::-1])
                if (b[1] != '.'):
                    sol_current += b[1]
                if k < len(str) - 1:
                    path_current.append((a, b, k+1, S_current))
                    solve(str, k + 1, b[0], S_current, sol_current, path_current)

        a = (q, '.', '.') #aleg tranzitiile care au lambda in loc de litera curenta din cuvant si nu scot simbolul din varful stivei
        # print(D[a])
        if a in D:
            for b in D[a]:
                #print(b)
                S_current = S.copy()
                sol_current = sol
                path_current = path.copy()

                # del S_current[-1]
                if (b[2] != '.'):
                    S_current += list(b[2][::-1])
                if (b[1] != '.'):
                    sol_current += b[1]
                path_current.append((a, b, k, S_current))
                solve(str, k, b[0], S_current, sol_current, path_current)


while True: #citesc cuvintele
    str = f.readline()[:-1]
    if not str:
        break
    print("Rezultatele pentru cuvantul " + str + ":")
    str += '.' #adaug . ca sa nu verific in algoritm daca am terminat cuvantul
    S = [] #stiva
    S.append(z0) #aduga simbolul de start al stivei
    sol = "" #cuvant output
    path = [] #drumul parcurs pana la starea finala (daca se ajunge)
    solve(str, 0, q0, S, sol, path) #va afisa toate cuvintele rezultate