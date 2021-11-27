from collections import deque
import sys
import copy

def add_dot(production,sign):
    left, right = production
    return [left, '.' + right + ';' + sign]

def move_dot(production):
    left, right = production
    pos = right.find('.')
    characters = list(right)
    characters[pos], characters[pos+1] = characters[pos+1], characters[pos]
    right = "".join(characters)
    return [left, right]

def next_chr_dot(production):
    left, right = production
    for i in range(len(right)):
        if right[i] == '.' and i < len(right) - 1:
            return right[i + 1]
    return ''

def make_first(secventa,grammar):
    res = []
    res.append("#")

    for i in range(len(secventa)):
        if secventa[i].isupper():
            for prod in grammar:
                left, right = prod
                if left == secventa[0]:
                    resY = make_first(right, grammar)
                    res.remove("#")
                    for rs in resY:
                        if rs not in res:
                            res.append(rs)
        elif secventa[i] not in res:
            if secventa[i] == "@":
                if "#" not in res:
                    res.append("#")
            else:
                res.append(secventa[i])

        if "#" not in res:
            return res

    return res

def first(production,grammar):
    # caut .neterminal
    # dau o lista cu toate terminalele
    left, right = production
    terminals = []
    for i in range(len(right)):
        if right[i] == '.' and i + 1 < len(right):
            if right[i+1].isupper() and i + 2 < len(right):
                quest = []
                for j in range(i+2, len(right)):
                    if right[j] != ";":
                        quest += right[j]
                return make_first(quest, grammar)

    return terminals

def closure(config, grammar):
    nonterminals = []
    for left, right in grammar:
        for c in left:
            if c.isupper() and c not in nonterminals:
                nonterminals.append(c)
        for c in right:
            if c.isupper() and c not in nonterminals:
                nonterminals.append(c)
    sol = list(config)
    dq = deque(sol)
    while len(dq) > 0:
        nod = dq.popleft()
        s = next_chr_dot(nod) # luam caracterul de dupa PUNCT
        if s in nonterminals: # daca acesta este nonterminal
            for prod in grammar: # ne uitam in productii
                if prod[0] == s: # daca productia incepe cu s
                    # pentru toate caracterele care sunt terminale
                    for i in first(nod, grammar):
                        nxt = add_dot(prod, i)
                        if nxt not in sol:
                            sol.append(nxt)
                            dq.append(nxt)
    return sol

def goto(config,nextChr,grammar):
    res = []
    for x in config:
        if nextChr == next_chr_dot(x):
            res.append(move_dot(x))
    return closure(res, grammar)

def config(grammar):
    symbols = []

    for left, right in grammar:
        for c in left:
            if c not in symbols:
                symbols.append(c)
        for c in right:
            if c not in symbols:
                symbols.append(c)

    configurations = [closure([add_dot(grammar[0], "#")], grammar)]
    transitions = {}
    start = 0
    sfarsit = 1
    while start < sfarsit:
        transitions[start] = [None] * len(symbols)
        for j in range(len(symbols)):
            nextConfig = goto(configurations[start], symbols[j], grammar)

            if not nextConfig:
                continue

            if nextConfig not in configurations:
                configurations.append(nextConfig)
                sfarsit += 1

            transitions[start][j] = configurations.index(nextConfig)

        start += 1

    return configurations, symbols, transitions

def eqq(a,b):
    aa, bb = a
    cc, dd = b
    if aa != cc:
        return False
    if bb != dd:
        return False
    return True

def noteqq(a,b):
    return not eqq(a,b)

def create_parsing_table(grammar):
    # calculez multimile canonice
    configurations, symbols, transitions = config(grammar)
    for i in range(len(configurations)):
        print("L" + str(i))
        for y in configurations[i]:
            print(y)
    symbols = []
    for left, right in grammar:
        for c in left:
            if c not in symbols and not c.isupper():
                symbols.append(c)
        for c in right:
            if c not in symbols and not c.isupper():
                symbols.append(c)

    if "#" not in symbols:
        symbols.append("#")

    for left, right in grammar:
        for c in left:
            if c not in symbols:
                symbols.append(c)
        for c in right:
            if c not in symbols:
                symbols.append(c)

    tabela = [dict() for i in configurations]
    for j in range(len(configurations)):
        for a in symbols:
            tabela[j][a] = ["error", 0]

    for j in range(len(configurations)):
        Ij = configurations[j]
        for prod in Ij:
            a = next_chr_dot(prod)
            if a != ";":

                Ik = goto(Ij, a, grammar)
                if Ik in configurations:
                    k = configurations.index(Ik)

                if a in symbols:
                    if a.isupper():
                        tabela[j][a] = ["goto",k]
                    else:
                        if eqq(tabela[j][a],["error", 0]):
                            tabela[j][a] = ["shift",k]
                        elif noteqq(tabela[j][a],["shift",k]):
                            print("Gramatica nu este LR(1)")
                            exit(0)
            else:
                left, right = prod
                if left == "Z" and right == "S.;#":
                    if eqq(tabela[j]["#"],["error", 0]):
                        tabela[j]["#"] = ["accept", 0]
                    elif noteqq(tabela[j]["#"] , ["accept", 0]):
                        print("Gramatica nu este LR(1)")
                        exit(0)
                else:
                    if eqq(tabela[j][right[-1]],["error", 0]):
                        tabela[j][right[-1]] = ["reduce", grammar.index([left,right[:-3]])]
                    elif noteqq(tabela[j][right[-1]],["reduce", grammar.index([left,right[:-3]])]):
                        print("Gramatica nu este LR(1)")
                        exit(0)

    print("Gramatica este LR(1)")
    return tabela

def LR_parser(tabela,grammar,query):
    characters = list(query)
    query += '#'
    characters.append('#')
    state_stack = [0]
    symbol_stack = []
    process = ""
    reductions = []
    index = 0

    actions = []
    actions.append("(0, " + query + ", " + "lambda)")
    not_done_yet = False
    while not_done_yet == False and index < len(characters):
        symbol = characters[index]
        s = state_stack[-1]

        if tabela[s][symbol] == ["error",0]:
            print("Eroare in timpul parsarii " + str(index) + " avand simbolul " + symbol)
            sys.exit()
        (action, value) = tabela[s][symbol]

        if action == "shift":
            symbol_stack.append(symbol)
            state_stack.append(value)
            index += 1
            nsymbol = copy.deepcopy(symbol)
            if nsymbol == '':
                nsymbol = 'lambda'
            process = "(Shiftare" + str(value) + " " + nsymbol + ")"
        elif action == "reduce":
            reductions.append(str(value))
            (left, right) = grammar[value]
            for i in range(len(right)):
                state_stack.pop()
                symbol_stack.pop()
            s = state_stack[-1]
            state_stack.append(tabela[s][left][1])
            symbol_stack.append(left)
            pleft = copy.deepcopy(left)
            pright = copy.deepcopy(right)
            if pleft == '':
                pleft = 'lambda'
            if pright == '':
                pright = 'lambda'
            process = "(Reducere" + str(value) + " " + pleft + " -> " + pright + ")"
        elif action == "accept":
            print("Query-ul {} este acceptat".format(query))
            for msg in actions:
                print(msg, end=" ->\n")
            print("accept")
            not_done_yet = True
        res = "(0"
        for i in range(len(symbol_stack)):
            res += str(symbol_stack[i]) + str(state_stack[i+1])
        res += ", " + query[index:]
        res += ", " + ''.join(reductions[::-1])
        if len(reductions) == 0:
            res += "lambda"
        res += ")"
        actions.append(process + " " + res)

def show_table(table,grammar):
    symbols = []
    for left, right in grammar:
        for c in left:
            if not c.isupper() and c not in symbols:
                symbols.append(c)
        for c in right:
            if not c.isupper() and c not in symbols:
                symbols.append(c)
    symbols.append('#')
    for left, right in grammar:
        for c in left:
            if c != 'Z' and c not in symbols:
                symbols.append(c)
        for c in right:
            if c != 'Z' and c not in symbols:
                symbols.append(c)

    print(symbols)
    print(" M",end=" ")
    for s in symbols:
        print("|   ",end=" ")
        print(s,end="    ")
    print()
    for k in range(len(table)):
        print(k,end=" " * (3 - len(str(k))))
        for s in symbols:
            print("|", end=" ")
            x = table[k][s][0] + str(table[k][s][1])
            print(x,end=" " * (8 - len(x)))
        print()

# lambda se regaseste ca find @ in input
def input_reader(inputFileName):
    f = open(inputFileName, "r")
    productions = []
    start_state = (f.readline())[:-1]
    productions.append(["Z", start_state])
    productions_number = int(f.readline()[:-1])
    for i in range(productions_number):
        x = (f.readline()[:-1]).split()
        left_symbol = x[0]
        right_symbol = x[2]
        if right_symbol == "@":
           right_symbol = ''
        productions.append([left_symbol, right_symbol])

    queries_number = int(f.readline()[:-1])
    queries = []
    for i in range(queries_number):
        query = (f.readline())[:-1]
        queries.append(query)

    return productions, queries

if __name__ == "__main__":
    grammar, queries = input_reader('input.txt')
    table = create_parsing_table(grammar)
    show_table(table,grammar)

    for query in queries:
        LR_parser(table,grammar,query)