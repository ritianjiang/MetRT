
import os
os.chdir('/home/owht/Download/data_from_web/temp')
gtf = open("Homo_sapiens.GRCh38.85.gtf")
gtfList = gtf.readlines()
aaa = file("LncRNA_merely_in_C_D1.result.txt","w+") 
count = 1

for line in open("LncRNA_merely_in_C_D1.csv"):
    print(count)
    infoList = line.split("\t")
    transName = infoList[4]
    for num in xrange(5,len(gtfList)):
        splitList = gtfList[num].split("\"")
        if splitList[5] == transName:
            
            nBack = num
            while ((gtfList[nBack].split("\""))[9]!='protein_coding'):
                if nBack > len(gtfList)-1:
                    break
                nBack+=1
            nForward = num
            while ((gtfList[nForward].split("\""))[9]!='protein_coding'):
                if nForward < 6:
                    break
                nForward-=1


            indx = int(infoList[1])
            indx2 = int(infoList[2])
            forward = gtfList[nForward].split("\t")
            forward = int(forward[4])
            back = gtfList[nBack].split("\t")
            back = int(back[3])
            if ((indx - forward) > (back - indx2)):
                temp = gtfList[nBack].split("\"")
                temp = temp[5]
                aaa.write(transName)
                aaa.write(";")
                aaa.write(temp)
                aaa.write("\n")
            else:
                temp = gtfList[nForward].split("\"")
                temp = temp[5]
                aaa.write(transName)
                aaa.write(";")
                aaa.write(temp)
                aaa.write("\n")
            break
    count +=1

aaa.close()



    
