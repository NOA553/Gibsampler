#Noa Babchick
#322985243
#the script gets as an input sequence reads, identifies overlaps between them, find the right order of them and then reconstruct the genomic sequence from the overlapping reads.
#python 3
import re
from math import *
from operator import itemgetter
import random

filename="DNA.txt"

def GetDNA(filename):
    list8=[]
    seq=open(filename,'r')
    for line in seq:
        list8.append("".join(line.split("\n")))
    return list8
dna=GetDNA(filename)
k=4
def FindMotifs(dna, k):
    list1=[]
    for i in dna:
        if (i!=""):

            num=random.randint(0,int(len(i)-k))
            list1.append((i[num:num+k]).upper())
            str1=""
    return list1
def CountOccurances(motifs):
    list2=[]
    for i in range(0,k):
        sumA=0
        sumT=0
        sumC=0
        sumG=0
        list3=[]
        for element in motifs:
            if(element[i]=="A"):
                sumA+=1
            if(element[i]=="T"):
                sumT+=1
            if(element[i]=="C"):
                sumC+=1
            if(element[i]=="G"):
                sumG+=1
        list3.append(sumA)
        list3.append(sumT)
        list3.append(sumC)
        list3.append(sumG)
        list2.append(list3)
    return list2



def CalcScore(motifs): 
    mat=CountOccurances(motifs)
    sum1=0
    for i in mat:
        max_value = max(i)
        sum2 = sum(i)
        sum1+=sum2-max_value #the number of times that the common nucleutide didnt apppear in its position
    return sum1

            
def FindProfile(motifs):
    mat=CountOccurances(motifs)
    n=4+(sum(mat[0]))
    for i in range(0,len(mat)):
        for j in range(0,4):
            mat[i][j]+=1
            mat[i][j]= mat[i][j]/n

    return mat
    
def ProfileProb(sequence, profile):
    list5=[]
    list4=[]
    prob=1
    for i in range(0,len(sequence)-k+1):
        list4.append(sequence[i:i+k])
    for element in list4:
        prob=1
        for j in range(0,k):
            if(element[j]=="A"):
               prob*=profile[j][0]
            if(element[j]=="T"):
               prob*=profile[j][1]
            if(element[j]=="C"):
               prob*=profile[j][2]
            if(element[j]=="G"):
               prob*=profile[j][3]
        list5.append(prob)
    return list5


def GibbsSampler(dna, k, N):
    
    motifs=FindMotifs(dna,k)
    bestMotifs=motifs
    bestScore=CalcScore(motifs)
    for i in range(1,N):
        num=random.randint(0,len(dna)-2)
        motifs2=bestMotifs
        del motifs2[num]
        profile=FindProfile(motifs2)
        sequence=dna[num]
        list7=ProfileProb(sequence, profile)
        maxim=list7[0]
        index=0
        for m in range(1,len(list7)):
            if(maxim<list7[m]):
                maxim=list7[m]
                index=m
        motifs2.insert(num,sequence[index:index+k])
        score=CalcScore(motifs2)
        if(bestScore>score):#if the bestScore is higgher than the score of the new motif list sign the new list as the best one
            bestScore=score
            bestMotifs=motifs2
        return bestMotifs,bestScore



def repeatGibbsSampler(dna, k, N, repeats):
    """ This function repeats the whole process in order to reach convergence. """
    bestMotifs = FindMotifs(dna,k)    
    bestScore = CalcScore(bestMotifs)
    
    for i in range(repeats):
        (motifs, score) = GibbsSampler(dna, k, N)
        if score < bestScore:
            bestMotifs = motifs
            bestScore = score
    return bestMotifs

print(repeatGibbsSampler(dna, k,10000000000,1000))
