import numpy as np
import random
N = 0 #numarul de cazuri posibile

def paralele(p1, p2, p3, p4, p5, p6):
    global N
    v1 = p3 - p1 #vectorul director pentru dreapta care contine punctele p1 si p3
    v2 = p2 - p1 #vectorul director pentru dreapta care contine punctele p1 si p2

    v3 = p6 - p4 #vectorul director pentru dreapta care contine punctele p4 si p6
    v4 = p5 - p4 #vectorul director pentru dreapta care contine punctele p4 si p5

    # produsul vectorial al 2 vectori cu aceeasi origine va returna vectorul director al planului format de dreptele corespunzatoare
    cp1 = np.cross(v1, v2)
    a1, b1, c1 = cp1 #coeficienti ai planului
    d1 = np.dot(cp1, p3) #calculam termenul liber punand unul din punctele din plan in ecuatia data de coeficientii de mai sus

    #analog pentru celalalt plan
    cp2 = np.cross(v3, v4)
    a2, b2, c2 = cp2
    d2 = np.dot(cp2, p6)

    #ecuatiile celor 2 plane
    A = np.array([a1, b1, c1, d1])
    B = np.array([a2, b2, c2, d2])

    #coeficientii unei ecuatii nu pot fi toti 0 pentru ca atunci nu am mai avea un plan
    if np.count_nonzero(A == 0) == 4 or np.count_nonzero(B == 0) == 4:
        return False

    N += 1
    #avem nevoie ca dreapta data de intersectia celor 2 plane sa fie paralela cu una din axe
    #deci orice puncte de pe aceasta dreapta trebuie sa aiba 2 coordonate constante
    #sa luam cazul in care dreapta este paralela cu Oz, deci coordonatele x si y sunt constante
    #pentru acest lucru vom considera ca z este necunoscuta secundara in sistemul dat de cele 2 ecuatii si o vom nota cu t
    #folosind metoda Cramer vom obtine ca pentru ca x si y sa fie constante (sa nu depinda de t)
    #avem nevoie de urmatoarele 2 ecuatii b1*c2*t = b2*c1*t si a1*c2*t = a2*c1*t pntru orice t
    #in cazul in care nu avem coeficienti nuli care sa satisfaca ecuatiile, va rezulta ca
    #a1/a2 = b1/b2 = c1/c2 lucru care inseamna ca cele 2 plane sunt paralele sau se suprapun, ceea ce nu ne intereseaza
    #pentru variantele in care nu avem c1 si c2 nule vom obtine pentru x sau y o valoare care depinde de t
    #din variantele in care c1 si c2 sunt nule le excludem pe cele in care a1 si a2 sau b1 si b2 sunt nule pentru ca nu am avea intersectie
    #asadar avem nevoie ca in matricea corespunzatoare sistemului coloana cu c sa fie 0 si niciuna din celelalte 2 sa fie 0
    #analog pentru cazurile in care drepta ar fi paralela cu alta axa
    #asadar avem nevoie ca in matricea sistemului sa existe o singura coloana 0
    k=0 #numar de coloane cu 0
    for i in range(3):
        if A[i]==B[i] and B[i]==0:
            k+=1
    if k==1:
        return True
    return False

s = 0 #numar cazuri favorabile
m=1 #numar cifre dupa virgula
l=1 #lungimea laturii cubului
n=100000 #numar incercari
#am considerat un cub cu un varf in origine si laturi paralele cu axele, orientate in sens pozitiv si de lungime l
#probabilitatea ramane aceeasi pentru orice cub de latura l
for i in range(n):
    p1 = np.array([round(random.random()*l, m), round(random.random()*l, m), round(random.random()*l, m)])
    p2 = np.array([round(random.random()*l, m), round(random.random()*l, m), round(random.random()*l, m)])
    p3 = np.array([round(random.random()*l, m), round(random.random()*l, m), round(random.random()*l, m)])
    p4 = np.array([round(random.random()*l, m), round(random.random()*l, m), round(random.random()*l, m)])
    p5 = np.array([round(random.random()*l, m), round(random.random()*l, m), round(random.random()*l, m)])
    p6 = np.array([round(random.random()*l, m), round(random.random()*l, m), round(random.random()*l, m)])
    if(paralele(p1,p2,p3,p4,p5,p6)):
        s += 1

print(s/N)
