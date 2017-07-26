import csv
import os
uuid_aliquote = csv.reader(open('/home/owht/KIZ/data/Deconvolution/XCI/aliquote_KIRC_uuid_bam_all.csv'))
uuid = list()
aliquote = list()
#这个csvreader对象，如果对它的row做for循环只能循环一次 - -、
for row in uuid_aliquote:
    uuid.append(row[1])
    aliquote.append(row[2])

#生成词典，每个uuid对应着一个aliquote的barcode号
dict_u_a = dict(zip(uuid,aliquote))

#设置当前目录为KIRC_with_uuid
os.chdir("KIRC_with_uuid")

filelist = os.listdir('.')
lenFL = len(filelist) - 1

i = 0
while i < lenFL:
    try:
        raw_name = filelist[i]
        file_uuid = raw_name.split(".")
        uid = file_uuid[1]
        print dict_u_a[uid]
        barcode = dict_u_a[uid]
        newname = raw_name.replace(uid,barcode)
        print newname
        os.rename(raw_name,newname) 
    except KeyError:
        print("the uuid is wrong!!")

    i+=1
